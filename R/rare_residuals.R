#' Detect Rare Events in Time Series Residuals
#'
#' This function filters a time series to compute residuals using STL decomposition
#' (for seasonal data) or LOESS smoothing (for non-seasonal data) and identifies
#' rare events (anomalies) in the residuals using Isolation Forest or DBSCAN. For
#' seasonal time series, it applies STL decomposition and full-length Fourier
#' smoothing separately to model periodic structure and estimate residuals.
#'
#' @param x A numeric vector (time series values), a univariate \code{ts} object, or
#'   a data frame with columns \code{time} (numeric or Date) and \code{value} (numeric).
#'   If a vector, assumes regular time steps starting from 1. If a \code{ts} is
#'   provided, the time index is derived via \code{stats::time(x)}.
#' @param seasonal Logical, indicating if the time series is seasonal (TRUE) or
#'   non-seasonal (FALSE) (default: FALSE). If \code{x} is a univariate \code{ts} object with
#'   \code{frequency(x) > 1}, the function will automatically treat the series as
#'   seasonal (i.e., set \code{seasonal = TRUE}).
#' @param period Numeric, the period of seasonality, e.g., 12 for monthly or 365.25
#'   for daily data (default: NULL). If NULL and \code{seasonal == TRUE}, the period is
#'   automatically determined. For \code{ts} input, \code{frequency(x)} is used; otherwise, it is
#'   estimated using spectral analysis.
#' @param method Character, specifying the anomaly detection method: "iforest",
#'   "dbscan", or "all" (default: "all").
#' @param fourier_terms Integer, number of Fourier term pairs for seasonal modeling
#'   in Fourier smoothing (default: 2).
#' @param stl_args A named list of arguments passed to \code{stats::stl()} for detrending and
#'    deseasonalizing \code{x}. Default: empty list \code{list()}.
#'    If not provided or empty, dynamic defaults for stable seasonal decomposition are used:
#'    \code{s.window} is the nearest odd to  \code{max(7, 1.5 * period)}, \code{t.window} is the nearest odd to
#'    \code{max(13, 1.5 * period)}, and \code{robust = TRUE}. If you provide a partial list,
#'    missing keys are filled with these dynamic defaults.
#' @param loess_args A named list of arguments passed to \code{stats::loess()} for
#'   non-seasonal smoothing. Default: \code{list()}. If not provided or empty,
#'   defaults \code{span = 0.75} and \code{degree = 1} are used. If you provide a
#'   partial list, missing keys are filled with these defaults. The formula is
#'   set to \code{value ~ time} unless you explicitly provide \code{formula}.
#' @param climatology_period Length-2 character vector with start and end dates
#'   (inclusive) for the climatology period used to estimate Fourier decomposition
#'   in seasonal data, e.g., c("1981-01-01","2010-12-31"). Only used when
#'   \code{seasonal = TRUE} and \code{time} is Date. Defaults to a 30-year period
#'   if available in the data; otherwise the full range of available dates is used.
#' @param iforest_args A named list of arguments passed to \code{rare_iforest()} (e.g.,
#'   \code{list(ntrees = 200)}). Default: \code{list()}.
#' @param dbscan_args A named list of arguments passed to \code{rare_dbscan()} (e.g.,
#'   \code{list(minPts = 10)}). Default: \code{list()}.
#' @return A list containing:
#'   \itemize{
#'     \item \code{data}: A data frame with:
#'       \itemize{
#'         \item \code{time}: Input time points.
#'         \item \code{value}: Original time series values.
#'         \item \code{residual}: Residuals from STL (seasonal) or LOESS (non-seasonal) smoothing.
#'         \item \code{residual_fourier}: Residuals from Fourier smoothing (if \code{seasonal = TRUE}).
#'         \item \code{is_anomaly_iforest}: Logical, anomalies from Isolation Forest (if applicable).
#'         \item \code{is_anomaly_dbscan}: Logical, anomalies from DBSCAN (if applicable).
#'         \item \code{score_iforest}: Anomaly scores from Isolation Forest.
#'         \item \code{score_dbscan}: Anomaly scores from DBSCAN.
#'       }
#'     \item \code{params}: A list of parameters used, including \code{period},
#'       \code{climatology_period} (if seasonal data with Date time index), and
#'       \code{method}.
#'   }
#' @details
#' For non-seasonal data, residuals are computed using LOESS with \code{time} as the
#' predictor. For seasonal data, residuals are computed twice: (1) using STL
#' decomposition to extract the remainder component, and (2) using a full-length
#' Fourier series to capture fixed periodicity. The STL seasonal component is
#' smoothed using the \code{s.window} parameter (numeric for flexible smoothing, "periodic"
#' for fixed seasonality). Rare events are detected in STL (or LOESS for non-seasonal)
#' residuals using \code{rare_iforest} or \code{rare_dbscan}. The period is estimated via
#' spectral analysis if not provided. Input validation prevents coercion errors.
#' Separate argument lists (\code{stl_args}, \code{iforest_args}, \code{dbscan_args}) ensure
#' function-specific parameters are passed correctly.
#'
#' For seasonal data with Date-based time indices, the \code{climatology_period} parameter
#' allows specification of a reference period for estimating the Fourier seasonal cycle.
#' The Fourier model is fitted using only data within this period. This is useful for
#' detecting changes in seasonality patterns relative to a baseline climatology.
#'
#' If \code{x} is a univariate \code{ts} object with \code{frequency(x) > 1}, the function will
#' automatically treat the series as seasonal (i.e., set \code{seasonal = TRUE}). When
#' \code{period} is not provided, \code{frequency(x)} will be used.
#'
#' @examples
#' \donttest{
#' # Annapolis temperature (seasonal, Date index)
#' data(Annapolis)
#' x <- data.frame(time = Annapolis$date, value = Annapolis$tmax)
#' result_tmax <- rare_residuals(x,
#'                          seasonal = TRUE, period = 365.25,
#'                          method = "iforest",
#'                          iforest_args = list(ntrees = 100))
#' plot(result_tmax$data$time, result_tmax$data$value, type = "l", las = 1,
#'      xlab = "Date", ylab = "Temperature (°C)",
#'      main = "Anomalies in Annapolis daily maximum temperature",
#' )
#' anomaly_sign <- sign(result_tmax$data$residual[result_tmax$data$is_anomaly_iforest])
#' points(result_tmax$data$time[result_tmax$data$is_anomaly_iforest],
#'        result_tmax$data$value[result_tmax$data$is_anomaly_iforest],
#'        col = c("blue", NA, "red")[anomaly_sign + 2],
#'        pch = 19)
#'
#' # Annapolis annual fish catch example (hypothetical data)
#' years <- 1990:2024
#' # Annual fish catch (tonnes)
#' fish_catch <- c(505, 493, 498, 487, 450, 479, 488, 472, 465, 476,
#'                 458, 466, 453, 461, 500, 447, 441, 435, 428, 443,
#'                 425, 418, 432, 421, 345, 412, 407, 418, 405, 397,
#'                 411, 403, 395, 408, 391)
#' result_fish <- rare_residuals(data.frame(time = years, value = fish_catch),
#'                               method = "iforest",
#'                               iforest_args = list(ntrees = 100))
#'
#' # Aggregate Annapolis daily data to identify rare warm spring years
#' # and test association with rare fish-catch years
#' spring_months <- c(3, 4, 5)
#' rare_tmax_spring <- result_tmax$data |>
#'     dplyr::mutate(Year = format(time, "%Y"),
#'                   Month = as.numeric(format(time, "%m"))) |>
#'     dplyr::filter(Month %in% spring_months, Year %in% years) |>
#'     dplyr::group_by(Year) |>
#'     dplyr::summarise(rare_tmax_spring = any(is_anomaly_iforest & residual > 0)) |>
#'     dplyr::arrange(Year) |>
#'     dplyr::pull(rare_tmax_spring)
#'
#' # Test association between rare fish-catch years and rare warm years
#' mcc_result <- mcc(result_fish$data$is_anomaly_iforest,
#'                   rare_tmax_spring,
#'                   ts = TRUE, bootstrap_reps = 999)
#' mcc_result$mcc
#' mcc_result$mcc_bootstrap_pv
#'
#' # Classic data example (ts input)
#' result_AirPassengers <- rare_residuals(AirPassengers,
#'                                       method = "iforest",
#'                                       iforest_args = list(ntrees = 100,
#'                                                           threshold = 0.6))
#' # View results
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot2::ggplot(result_AirPassengers$data, aes(x = time, y = value)) +
#'      geom_line(color = "gray") +
#'      geom_point(aes(color = is_anomaly_iforest)) +
#'      labs(title = "Anomaly detection in AirPassengers using isolation forest",
#'           x = "Time",
#'           y = "Number of passengers",
#'           color = "Anomaly status") +
#'      scale_color_manual(values = c("gray", "red"), labels = c("Normal", "Anomaly")) +
#'      theme_minimal()
#' }
#' }
#' @seealso \code{\link{rare_iforest}}, \code{\link{rare_dbscan}},  \code{\link{mcc}}, \code{\link[stats]{stl}}
#' @importFrom stats loess na.pass predict spec.pgram lm stl ts
#' @export
#'
rare_residuals <- function(x, seasonal = FALSE, period = NULL,
                           method = c("all", "iforest", "dbscan"),
                           fourier_terms = 2,
                           stl_args = list(),
                           loess_args = list(),
                           climatology_period = NULL,
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

    # If input is ts with frequency > 1, auto-enable seasonality (unless already TRUE)
    if (!is.null(original_ts)) {
        ts_freq <- stats::frequency(original_ts)
        if (is.finite(ts_freq) && ts_freq > 1 && !isTRUE(seasonal)) {
            seasonal <- TRUE
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

    # Determine climatology period for Fourier decomposition (seasonal data only)
    clim_indices <- NULL
    if (seasonal && inherits(time, "Date")) {
        # Default climatology period: prefer 30-year window if available, else full range
        d_min <- min(time, na.rm = TRUE)
        d_max <- max(time, na.rm = TRUE)
        if (is.null(climatology_period)) {
            start_year <- as.integer(format(d_min, "%Y"))
            end_year <- as.integer(format(d_max, "%Y"))
            if (end_year - start_year + 1 >= 30) {
                climatology_period <- c(sprintf("%d-01-01", start_year), sprintf("%d-12-31", start_year + 29))
            } else {
                climatology_period <- c(as.character(d_min), as.character(d_max))
            }
        }
        # Find indices within climatology period
        clim_start <- as.Date(climatology_period[1])
        clim_end <- as.Date(climatology_period[2])
        clim_indices <- which(time >= clim_start & time <= clim_end)
        if (length(clim_indices) == 0) {
            warning("No data found in specified climatology_period. Using full time series for Fourier decomposition.")
            clim_indices <- NULL
        }
    }

    # Smoothing and residual computation
    if (seasonal) {
        # STL decomposition for seasonal data
        ts_data <- stats::ts(value, frequency = period)
        # Dynamic defaults for stable seasonal decomposition
        nearest_odd <- function(z, min_val) {
            z <- as.integer(round(z))
            if (z < min_val) z <- min_val
            if (z %% 2 == 0) z <- z + 1L
            z
        }
        period_int <- max(2L, as.integer(round(period)))
        # Propose windows based on period and length
        s_prop <- nearest_odd(max(7L, as.integer(round(1.5 * period_int))), 7L)
        t_prop <- nearest_odd(max(13L, as.integer(round(1.5 * period_int))), 3L)
        # Clamp to series length constraints
        s_prop <- min(s_prop, if (n %% 2 == 0) n - 1L else n - 2L)
        if (s_prop < 7L) s_prop <- 7L
        t_prop <- min(t_prop, if (n %% 2 == 0) n - 1L else n - 2L)
        if (t_prop < 3L) t_prop <- 3L

        stl_args_mod <- stl_args
        # If user omitted stl_args or provided an empty list, apply dynamic defaults
        if (missing(stl_args) || is.null(stl_args_mod) || length(stl_args_mod) == 0) {
            stl_args_mod <- list(s.window = s_prop, t.window = t_prop, robust = TRUE)
        } else {
            if (is.null(stl_args_mod$s.window)) stl_args_mod$s.window <- s_prop
            if (is.null(stl_args_mod$t.window)) stl_args_mod$t.window <- t_prop
            if (is.null(stl_args_mod$robust))   stl_args_mod$robust   <- TRUE
        }

        stl_fit <- do.call(stats::stl, c(list(x = ts_data), stl_args_mod))
        stl_seasonal <- stl_fit$time.series[, "seasonal"]
        # Additional smoothing of the seasonal cycle estimated on daily data
        # if (period > 12) {
        #     # Make span proportional to one year of the data, period / n, with limits
        #     span <- max(0.1, min(0.9, period / n))
        #     stl_seasonal <- stats::loess(stl_seasonal ~ seq_along(stl_seasonal), span = span)$fitted
        # }
        # Reconstruct fitted values and residuals
        stl_fitted <- stl_seasonal + stl_fit$time.series[, "trend"]
        residual <- ts_data - stl_fitted

        # Compute Fourier terms
        i <- 0:(n - 1)
        omega <- 2 * pi * i / period
        fourier <- data.frame(
            matrix(0, nrow = n, ncol = 2 * fourier_terms)
        )
        colnames(fourier) <- paste0(rep(c("sin", "cos"), each = fourier_terms), 1:fourier_terms)
        for (k in 1:fourier_terms) {
            fourier[[paste0("sin", k)]] <- sin(omega * k)
            fourier[[paste0("cos", k)]] <- cos(omega * k)
        }

        # Fourier smoothing: fit on climatology period (if specified), predict on full series
        if (!is.null(clim_indices)) {
            # Fit Fourier model on climatology period only
            fourier_data_clim <- data.frame(
                value = (value - stl_fit$time.series[, "trend"])[clim_indices],
                fourier[clim_indices, , drop = FALSE]
            )
            fourier_fit <- stats::lm(value ~ ., data = fourier_data_clim)
            # Predict on full time series
            fourier_data_full <- data.frame(fourier)
            fourier_fitted <- stats::predict(fourier_fit, newdata = fourier_data_full)
        } else {
            # Fit and predict on full time series
            fourier_data <- data.frame(value = value - stl_fit$time.series[, "trend"], fourier)
            fourier_fit <- stats::lm(value ~ ., data = fourier_data)
            fourier_fitted <- stats::predict(fourier_fit)
        }
        residual_fourier <- value - fourier_fitted
    } else {
        # LOESS for non-seasonal data
        loess_args_mod <- loess_args
        if (missing(loess_args) || is.null(loess_args_mod) || length(loess_args_mod) == 0) {
            loess_args_mod <- list(formula = value ~ time, span = 0.75, degree = 1)
        } else {
            if (is.null(loess_args_mod$formula)) loess_args_mod$formula <- value ~ time
            if (is.null(loess_args_mod$span))    loess_args_mod$span    <- 0.75
            if (is.null(loess_args_mod$degree))  loess_args_mod$degree  <- 1
        }
        loess_fit <- do.call(stats::loess, loess_args_mod)
        loess_fitted <- stats::predict(loess_fit)
        residual <- value - loess_fitted
        residual_fourier <- rep(NA, n)
    }

    # Prepare output data frame
    output <- data.frame(
        time = time,
        value = value
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
    params <- list(
        period = if (seasonal) period else NA,
        method = method
    )
    if (seasonal && inherits(time, "Date") && !is.null(climatology_period)) {
        params$climatology_period <- climatology_period
    }

    return(list(
        data = output,
        params = params
    ))
}
