#' Detect Rare Events and Seasonality Shifts in Time Series Residuals
#'
#' This function filters a time series to compute residuals using STL decomposition
#' (for seasonal data) or LOESS smoothing (for non-seasonal data) and identifies
#' rare events (anomalies) in the residuals using Isolation Forest or DBSCAN. For
#' seasonal time series, it applies STL decomposition and full-length Fourier
#' smoothing separately, comparing their seasonality estimates to detect shifts
#' (e.g., early spring) by analyzing differences in seasonal peaks.
#'
#' @param x A numeric vector (time series values), a univariate `ts` object, or
#'   a data frame with columns `time` (numeric or Date) and `value` (numeric).
#'   If a vector, assumes regular time steps starting from 1. If a `ts` is
#'   provided, the time index is derived via `stats::time(x)`.
#' @param seasonal Logical, indicating if the time series is seasonal (TRUE) or
#'   non-seasonal (FALSE) (default: FALSE). If `x` is a univariate `ts` object with
#'   `frequency(x) > 1`, the function will automatically treat the series as
#'   seasonal (i.e., set `seasonal = TRUE`).
#' @param period Numeric, the period of seasonality, e.g., 12 for monthly or 365.25
#'   for daily data (default: NULL). If NULL and `seasonal == TRUE`, the period is
#'   automatically determined. For `ts` input, `frequency(x)` is used; otherwise, it is
#'   estimated using spectral analysis.
#' @param method Character, specifying the anomaly detection method: "iforest",
#'   "dbscan", or "all" (default: "all").
#' @param fourier_terms Integer, number of Fourier term pairs for seasonal modeling
#'   in Fourier smoothing (default: 2).
#' @param seasonality_shift Logical, whether to compute seasonality shift (e.g., early
#'   spring) by comparing STL and Fourier seasonal components (default: TRUE if
#'   `seasonal = TRUE`).
#' @param stl_args A named list of arguments passed to `stats::stl()` for detrending and
#'    deseasonalizing `x`. Default: `list(s.window = 7, robust = TRUE)`.
#' @param iforest_args A named list of arguments passed to `rare_iforest()` (e.g.,
#'   `list(ntrees = 200)`). Default: `list()`.
#' @param dbscan_args A named list of arguments passed to `rare_dbscan()` (e.g.,
#'   `list(minPts = 10)`). Default: `list()`.
#' @return A list containing:
#'   \itemize{
#'     \item \code{data}: A data frame with:
#'       \itemize{
#'         \item \code{time}: Input time points.
#'         \item \code{value}: Original time series values.
#'         \item \code{year}: Year extracted from time.
#'         \item \code{residual}: Residuals from STL (seasonal) or LOESS (non-seasonal) smoothing.
#'         \item \code{residual_fourier}: Residuals from Fourier smoothing (if `seasonal = TRUE`).
#'         \item \code{is_anomaly_iforest}: Logical, anomalies from Isolation Forest (if applicable).
#'         \item \code{is_anomaly_dbscan}: Logical, anomalies from DBSCAN (if applicable).
#'         \item \code{score_iforest}: Anomaly scores from Isolation Forest.
#'         \item \code{score_dbscan}: Anomaly scores from DBSCAN.
#'       }
#'     \item \code{seasonality_shift}: A list (if `seasonality_shift = TRUE` and
#'       `seasonal = TRUE`) with:
#'       \itemize{
#'         \item \code{lag}: Estimated time shift (positive if STL peaks
#'           earlier than Fourier, e.g., early spring).
#'         \item \code{stl_seasonal}: Seasonal component from STL.
#'         \item \code{fourier_seasonal}: Seasonal component from Fourier.
#'       }
#'   }
#' @details
#' For non-seasonal data, residuals are computed using LOESS with `year` as the
#' predictor. For seasonal data, residuals are computed twice: (1) using STL
#' decomposition to extract the remainder component, and (2) using a full-length
#' Fourier series to capture fixed periodicity. The STL seasonal component is
#' smoothed using the `s.window` parameter (numeric for flexible smoothing, "periodic"
#' for fixed seasonality) and further smoothed with LOESS for shift detection.
#' Seasonality shifts are detected by comparing smoothed STL and Fourier seasonal
#' components via cross-correlation to estimate the lag (in days) where peaks differ
#' (e.g., early spring). Rare events are detected in STL (or LOESS for non-seasonal)
#' residuals using `rare_iforest` or `rare_dbscan`. The period is estimated via
#' spectral analysis if not provided. Input validation prevents coercion errors.
#' Separate argument lists (`stl_args`, `iforest_args`, `dbscan_args`) ensure
#' function-specific parameters are passed correctly.
#'
#' If `x` is a univariate `ts` object with `frequency(x) > 1`, the function will
#' automatically treat the series as seasonal (i.e., set `seasonal = TRUE`). When
#' `period` is not provided, `frequency(x)` will be used.
#'
#' @examples
#' \dontrun{
#' # Seasonal time series with Date index
#' set.seed(123)
#' dates <- seq(as.Date("2020-01-01"), by = "day", length.out = 1500)
#' value <- sin(2 * pi * seq_along(dates) / 365.25) + rnorm(1500, 0, 0.2)
#' x <- data.frame(time = dates, value = value)
#' result <- rare_residuals(x, seasonal = TRUE, period = 365.25, method = "all",
#' iforest_args = list(ntrees = 100))
#' plot(result$data$time, result$data$value, type = "l", main = "Time Series")
#' points(result$data$time[result$data$is_anomaly_iforest],
#'        result$data$value[result$data$is_anomaly_iforest], col = "red", pch = 19)
#' cat("Seasonality shift:", result$seasonality_shift$lag, "\n")
#'
#' # Plot seasonal components
#' plot(result$data$time, result$seasonality_shift$stl_seasonal, type = "l",
#'      col = "blue", main = "Seasonal Components")
#' lines(result$data$time, result$seasonality_shift$fourier_seasonal, col = "red")
#' legend("topright", c("STL", "Fourier"), col = c("blue", "red"), lty = 1)
#'
#' # Classic data example (ts input)
#'result_AirPassengers <- rare_residuals(AirPassengers,
#'                                       method = "iforest",
#'                                       iforest_args = list(ntrees = 100,
#'                                                           threshold = 0.6))
#'# View results
#'ggplot2::ggplot(result_AirPassengers$data, aes(x = time, y = value)) +
#'    geom_line(color = "gray") +
#'    geom_point(aes(color = is_anomaly_iforest)) +
#'    labs(title = "Anomaly detection in AirPassengers using isolation forest",
#'         x = "Time",
#'         y = "Number of passengers",
#'         color = "Anomaly status") +
#'    scale_color_manual(values = c("gray", "red"), labels = c("Normal", "Anomaly")) +
#'    theme_minimal()
#' }
#' @importFrom stats loess na.pass predict spec.pgram lm ccf stl ts
#' @export
#'
rare_residuals <- function(x, seasonal = FALSE, period = NULL,
                           method = c("all", "iforest", "dbscan"),
                           fourier_terms = 2,
                           seasonality_shift = seasonal,
                           stl_args = list(s.window = 7, robust = TRUE),
                           iforest_args = list(),
                           dbscan_args = list()) {
    # Input validation
    original_ts <- NULL
    if (inherits(x, "ts")) {
        # Ensure univariate ts
        if (NCOL(as.matrix(x)) != 1) {
            stop("If `x` is a ts object, it must be univariate.")
        }
        original_ts <- x
        time <- stats::time(x)
        value <- as.numeric(x)
        x <- data.frame(time = time, value = value)
    } else if (is.vector(x) && is.numeric(x)) {
        time <- seq_along(x)
        value <- x
        x <- data.frame(time = time, value = value)
    } else if (is.data.frame(x)) {
        if (!all(c("time", "value") %in% colnames(x))) {
            stop("If `x` is a data frame, it must contain `time` and `value` columns.")
        }
        if (!is.numeric(x$value)) {
            stop("`value` column must be numeric.")
        }
        if (!(is.numeric(x$time) || inherits(x$time, "Date"))) {
            stop("`time` column must be numeric or Date.")
        }
    } else {
        stop("`x` must be a numeric vector, a univariate ts, or a data frame with `time` and `value` columns.")
    }

    if (any(!is.finite(x$value))) {
        stop("`value` contains NA, NaN, or Inf values. Please clean the data.")
    }
    # if (any(duplicated(x$time))) {
    #     stop("`time` contains duplicate values, which can cause numerical instability.")
    # }
    if (!is.logical(seasonal)) {
        stop("`seasonal` must be logical (TRUE or FALSE).")
    }
    method <- match.arg(method)
    if (seasonal && !is.null(period) && (!is.numeric(period) || period <= 0)) {
        stop("`period` must be a positive numeric value or NULL.")
    }
    if (!is.integer(fourier_terms) && !is.numeric(fourier_terms) || fourier_terms < 1) {
        stop("`fourier_terms` must be a positive integer.")
    }
    fourier_terms <- as.integer(fourier_terms)
    if (!is.logical(seasonality_shift)) {
        stop("`seasonality_shift` must be logical (TRUE or FALSE).")
    }

    # Check package availability
    if (method %in% c("iforest", "all") && !requireNamespace("isotree", quietly = TRUE)) {
        stop("Package 'isotree' is required for Isolation Forest. Install it using install.packages('isotree').")
    }
    if (method %in% c("dbscan", "all") && !requireNamespace("dbscan", quietly = TRUE)) {
        stop("Package 'dbscan' is required for DBSCAN. Install it using install.packages('dbscan').")
    }

    # Extract time and value
    time <- x$time
    value <- x$value
    n <- length(value)
    if (inherits(time, "Date")) {
        year <- as.numeric(format(time, "%Y"))
    } else {
        year <- time
    }

    # If input is ts with frequency > 1, auto-enable seasonality (unless already TRUE)
    if (!is.null(original_ts)) {
        ts_freq <- stats::frequency(original_ts)
        if (is.finite(ts_freq) && ts_freq > 1 && !isTRUE(seasonal)) {
            seasonal <- TRUE
            if (missing(seasonality_shift)) seasonality_shift <- TRUE
            message("Detected 'ts' with frequency > 1; setting seasonal = TRUE.")
        }
    }

    # Estimate/derive period if not provided for seasonal data
    if (seasonal && is.null(period)) {
        if (!is.null(original_ts)) {
            period <- stats::frequency(original_ts)
            message("`period` not specified. Using ts frequency: ", period, ".")
        } else {
            spec <- stats::spec.pgram(value, plot = FALSE)
            period <- 1 / spec$freq[which.max(spec$spec)]
            message("`period` not specified. Estimated as ", round(period, 2),
                    " using spectral analysis.")
        }
    }

    # Smoothing and residual computation
    if (seasonal) {
        # STL decomposition for seasonal data
        ts_data <- stats::ts(value, frequency = period)
        stl_fit <- do.call(stats::stl, c(list(x = ts_data), stl_args))
        stl_seasonal <- stl_fit$time.series[, "seasonal"]
        # Additional smoothing of the seasonal cycle estimated on daily data
        if (period > 12) {
            stl_seasonal <- stats::loess(stl_seasonal ~ seq_along(stl_seasonal), span = period/n)$fitted
        }
        # Reconstruct fitted values and residuals
        stl_fitted <- stl_seasonal + stl_fit$time.series[, "trend"]
        residual <- ts_data - stl_fitted

        # Compute Fourier terms
        t <- seq(0, 1, length.out = n) * (n - 1) / period
        fourier <- data.frame(
            matrix(0, nrow = n, ncol = 2 * fourier_terms)
        )
        colnames(fourier) <- paste0(rep(c("sin", "cos"), each = fourier_terms), 1:fourier_terms)
        for (k in 1:fourier_terms) {
            fourier[[paste0("sin", k)]] <- sin(2 * pi * k * t)
            fourier[[paste0("cos", k)]] <- cos(2 * pi * k * t)
        }

        # Fourier smoothing of detrended data
        fourier_data <- data.frame(value = value - stl_fit$time.series[, "trend"], fourier)
        fourier_fit <- stats::lm(value ~ ., data = fourier_data)
        fourier_fitted <- stats::predict(fourier_fit)
        residual_fourier <- value - fourier_fitted
        fourier_seasonal <- fourier_fitted
    } else {
        # LOESS for non-seasonal data
        loess_fit <- stats::loess(value ~ year, span = 0.75, degree = 1)
        loess_fitted <- stats::predict(loess_fit)
        residual <- value - loess_fitted
        residual_fourier <- rep(NA, n)
        fourier_fitted <- rep(NA, n)
        fourier_seasonal <- NULL
    }

    # Seasonality shift detection (for seasonal data)
    seasonality_shift_result <- NULL
    if (seasonal && seasonality_shift) {
        # Compute cross-correlation to estimate shift
        ccf_result <- stats::ccf(as.vector(stl_seasonal), fourier_seasonal, lag.max = round(period / 2),
                                 plot = FALSE, na.action = na.pass)
        lag <- ccf_result$lag[which.max(abs(ccf_result$acf))]
        seasonality_shift_result <- list(
            lag = lag,
            stl_seasonal = stl_seasonal,
            fourier_seasonal = fourier_seasonal
        )
    }

    # Prepare output data frame
    output <- data.frame(
        time = time,
        value = value,
        year = year
    )
    output$residual <- residual
    output$residual_fourier <- residual_fourier

    # Rare event detection
    residual_matrix <- matrix(residual, ncol = 1)

    # Isolation Forest
    if (method %in% c("iforest", "all")) {
        iforest_result <- do.call(rare_iforest, c(list(x = residual_matrix), iforest_args))
        output$is_anomaly_iforest <- iforest_result$is_anomaly
        output$score_iforest <- iforest_result$scores
    }

    # DBSCAN
    if (method %in% c("dbscan", "all")) {
        dbscan_result <- do.call(rare_dbscan, c(list(x = residual_matrix), dbscan_args))
        output$is_anomaly_dbscan <- dbscan_result$is_anomaly
        output$score_dbscan <- dbscan_result$scores
    }

    # Return results
    return(list(
        data = output,
        seasonality_shift = seasonality_shift_result
    ))
}
