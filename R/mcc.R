#' Calculate and Test the Matthews Correlation Coefficient (MCC)
#'
#' Based on the vectors of actual and predicted binary values, compute the confusion
#' matrix and the MCC (mean square contingency coefficient).
#' The function implements two methods for testing significance:
#' 1. A parametric chi-square test.
#' 2. A non-parametric bootstrap test for a more robust p-value.
#'
#' @param actual A vector of true or observed binary labels. Can be numeric (0/1),
#'   logical (TRUE/FALSE), character, or factor.
#' @param predicted A vector of predicted binary labels, of the same type and
#'   length as `actual`.
#' @param positive_class An optional value explicitly specifying the "positive"
#'   class label. If `NULL`, the function will infer it.
#' @param bootstrap_reps The number of bootstrap replicates for p-value calculation.
#'   Default is 999. Set to 0 to disable bootstrapping.
#' @param ts A logical flag indicating if the data are time series. Default is FALSE. If TRUE,
#'  a version of bootstrap for time series is applied to account for potential autocorrelation;
#'  the data are then assumed to be ordered in time.
#' @param confidence The confidence level for the bootstrap confidence interval (default: 0.95).
#'  Must be between 0 and 1 (exclusive).
#' @param l The block length for the time series bootstrap (default: 5). Only used if `ts` is TRUE,
#'  see `?boot::tsboot`.
#' @param sim The type of simulation for the time series bootstrap (default: "fixed").
#'  Options are "fixed" (moving block bootstrap) or "geom" (stationary bootstrap), see `?boot::tsboot`.
#' @param ... Additional arguments passed to `boot::tsboot` when `ts` is TRUE.
#'
#'
#' @return A list containing:
#'   - `confusion_matrix`: The 2x2 confusion matrix.
#'   - `mcc`: The Matthews Correlation Coefficient (-1 to +1).
#'   - `chi_square_test`: The output of `chisq.test()`. `$p.value` is the parametric p-value.
#'   - `mcc_bootstrap_pv`: The p-value from the bootstrap permutation test.
#'   - `mcc_bootstrap_ci`: The bootstrap confidence interval for the MCC.
#'   - `mcc_bootstrap_reps`: The number of bootstrap replicates used.
#'   - `positive_class`: The positive class label used.
#'   - `confidence`: The confidence level for the bootstrap confidence interval.
#'   - `ts`: Whether the data were treated as time series.
#'
#'
#' @details
#' The MCC is a robust metric for binary classification, especially on imbalanced data.
#' It ranges from -1 (perfect negative correlation) to +1 (perfect positive correlation).
#' For 2x2 confusion matrices:
#' \[
#' \text{MCC} = \frac{TP \cdot TN - FP \cdot FN}{\sqrt{(TP + FP)(TP + FN)(TN + FP)(TN + FN)}}
#' \]
#' where:
#' - TP: True Positives
#' - TN: True Negatives
#' - FP: False Positives
#' - FN: False Negatives
#'
#' and
#' \[
#' |\text{MCC}| = \sqrt{\frac{\chi^2}{n}}
#' \]
#'
#' The **chi-square test** is fast but assumes that each observation is independent,
#' which is often not true for time series data.
#'
#' **For Time Series:** Both standard chi-square and standard bootstrapping assume
#' data independence. If the `ts` parameter is set to TRUE, the function applies a
#' version of bootstrapping suitable for time series data to account for potential
#' autocorrelation. The data are assumed to be ordered in time. The block length `l`
#' can be adjusted based on the expected autocorrelation structure. The chi-square
#' test is still provided for reference but should be interpreted with caution.
#' The bootstrap p-value is more reliable for time series data.
#'
#' The **bootstrap test** is more computationally intensive but does not assume
#' independence and is more robust, especially for small sample sizes or when the
#' data may be autocorrelated.
#' The bootstrap p-value is calculated as the proportion of bootstrap MCC values
#' that are as extreme or more extreme than the observed MCC, using a two-tailed
#' test.
#' Confidence intervals for the MCC are also provided based on the bootstrap distribution.
#' The function uses the `boot` package for bootstrapping.
#'
#' @examples
#' # Example 1: A clear, significant correlation
#' actual_vals <- rep(c(1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1), 3)
#' predicted_vals <- rep(c(1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0), 3)
#' mcc_results <- mcc(actual_vals, predicted_vals, positive_class = 1)
#' print(mcc_results$confusion_matrix)
#' print(paste("MCC:", round(mcc_results$mcc, 3)))
#' print(paste("Chi-Square p-value:", round(mcc_results$chi_square_test$p.value, 4)))
#' print(paste("Bootstrap p-value:", round(mcc_results$mcc_bootstrap_pv, 4)))
#'
#' # Example 2: No significant correlation
#' actual_rand <- rep(c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0), 3)
#' predicted_rand <- rep(c(1, 1, 0, 0, 1, 1, 0, 0, 1, 0), 3)
#' mcc_rand <- mcc(actual_rand, predicted_rand, bootstrap_reps = 999)
#' print(paste("MCC:", round(mcc_rand$mcc, 3)))
#' print(paste("Bootstrap p-value:", round(mcc_rand$mcc_bootstrap_pv, 4)))
#'
#' # Example 3: Time series data with autocorrelation in both series
#' n <- 100
#' set.seed(12345)
#' time_series_actual <- rbinom(n, 1, 0.3)
#' time_series_predicted <- ifelse(runif(n) < 0.3,
#'                                time_series_actual,
#'                                1 - time_series_actual)
#' # Introduce autocorrelation
#' for (i in 2:n) {
#'    if (runif(1) < 0.7) {
#'        time_series_actual[i] <- time_series_actual[i - 1]
#'        time_series_predicted[i] <- time_series_predicted[i - 1]
#'    }
#' }
#'
#' # Check autocorrelation at lag 1
#' par(mfrow = c(1, 2))
#' mcc(time_series_actual, dplyr::lag(time_series_actual), ts = TRUE, l = 7)$mcc
#' acf(time_series_actual, main = "ACF of Actual Time Series, treated as numeric")
#' mcc(time_series_predicted, dplyr::lag(time_series_predicted), ts = TRUE, l = 7)$mcc
#' acf(time_series_predicted, main = "ACF of Predicted Time Series")
#'
#' # Visualize the time series
#' plot.ts(cbind(time_series_actual, time_series_predicted),
#'        main = "Time Series of Actual actual Predicted", xlab = "Time")
#'
#' # Calculate MCC treating data as time series
#' mcc_ts <- mcc(time_series_actual,
#'               time_series_predicted,
#'               positive_class = 1,
#'               ts = TRUE, l = 7, bootstrap_reps = 999)
#' print(paste("MCC:", round(mcc_ts$mcc, 3)))
#' print(paste("Chi-Square p-value:", round(mcc_ts$chi_square_test$p.value, 4)))
#' print(paste("Bootstrap p-value:", round(mcc_ts$mcc_bootstrap_pv, 4)))
#' print(paste("Bootstrap CI:",
#'    paste(round(mcc_ts$mcc_bootstrap_ci, 3), collapse = " to ")))
#'
#' @importFrom stats chisq.test quantile
#' @importFrom boot tsboot
#' @export
#'
mcc <- function(actual,
                predicted,
                positive_class = NULL,
                bootstrap_reps = 999,
                confidence = 0.95,
                ts = FALSE,
                l = 5,
                sim = "fixed", ...) {

  # --- Internal helper function to compute MCC from a matrix ---
  calculate_mcc_from_matrix <- function(mat) {
    TP <- as.numeric(mat[1, 1])
    FN <- as.numeric(mat[1, 2])
    FP <- as.numeric(mat[2, 1])
    TN <- as.numeric(mat[2, 2])
    numerator <- TP * TN - FP * FN
    denominator_parts <- c(TP + FP, TP + FN, TN + FP, TN + FN)
    if (any(denominator_parts == 0)) {
      # If any part of the denominator is zero, return 0 to avoid division by zero
      return(0)
    }
    denominator <- sqrt(prod(denominator_parts))
    numerator / denominator
  }

  # --- Internal helper function to compute MCC for bootstrapping ---
  calculate_mcc <- function(data, ordered_levels) {
    mat <- table(Actual = factor(c(data[, 1], ordered_levels, rev(ordered_levels)), levels = ordered_levels),
                 Predicted = factor(c(data[, 2], ordered_levels, ordered_levels), levels = ordered_levels)) - 1
    calculate_mcc_from_matrix(mat)
  }

  # --- Input Validation ---
  if (length(actual) != length(predicted)) {
    stop("Input vectors 'actual' and 'predicted' must have the same length.")
  }
  if (length(actual) == 0) {
    stop("Input vectors must not be empty.")
  }
  unique_values <- sort(unique(c(as.character(actual), as.character(predicted))))
  if (length(unique_values) < 2) {
      warning("Only one class is present in the data. MCC and related stats cannot be calculated.")
      # Return the standard output structure with NAs for calculated values
      return(
          list(
              confusion_matrix = NA,
              chi_square_test = NA,
              mcc = NA,
              mcc_bootstrap_pv = NA,
              mcc_bootstrap_ci = NA,
              mcc_bootstrap_reps = bootstrap_reps,
              positive_class = if (!is.null(positive_class)) positive_class else unique_values[1],
              confidence = confidence,
              ts = ts
          )
      )
  }
  if (length(unique_values) > 2) {
    warning(paste("Data contains more than two unique values (must be binary):",
                  paste(unique_values, collapse = ", ")))
  }
  if (length(unique(actual)) > 2) {
    stop("'actual' vector must be binary (contain at most two unique values).")
  }
  if (!is.logical(ts) || length(ts) != 1) {
    stop("`ts` must be a single logical value (TRUE or FALSE).")
  }
  if (!is.numeric(bootstrap_reps) || length(bootstrap_reps) != 1 ||
      bootstrap_reps < 0 || bootstrap_reps != round(bootstrap_reps)) {
    stop("`bootstrap_reps` must be a single non-negative integer.")
  }
  if (!is.numeric(confidence) || length(confidence) != 1 ||
      confidence <= 0 || confidence >= 1) {
    stop("`confidence` must be a single numeric value between 0 and 1 (0 < confidence < 1).")
  }
  if (bootstrap_reps > 0 && !requireNamespace("boot", quietly = TRUE)) {
    stop("Package 'boot' is required for bootstrapping. Please install it using install.packages('boot').")
  }

  # --- Determine Positive/Negative Classes ---
  if (is.null(positive_class)) {
    if (is.logical(actual)) {
      positive_class <- TRUE
    } else if (is.numeric(actual) && all(actual %in% c(0, 1))) {
      positive_class <- 1
    } else {
      positive_class <- unique_values[2]
    }
    warning(paste("No positive_class specified. Inferred positive class:",
                  positive_class))
  }
  positive_class <- as.character(positive_class)
  negative_class <- unique_values[unique_values != positive_class][1]

  # --- Create Confusion Matrix ---
  ordered_levels <- c(positive_class, negative_class)
  conf_matrix <- table(Actual = factor(c(actual, ordered_levels, rev(ordered_levels)), levels = ordered_levels),
                       Predicted = factor(c(predicted, ordered_levels, ordered_levels), levels = ordered_levels)) - 1

  if (any(rowSums(conf_matrix) == 0) || any(colSums(conf_matrix) == 0)) {
      warning("The confusion matrix has a row or column with all zeros. MCC is undefined and returned as NA.")
      return(
          list(
              confusion_matrix = conf_matrix, # Return the matrix for inspection
              chi_square_test = NA,
              mcc = NA,
              mcc_bootstrap_pv = NA,
              mcc_bootstrap_ci = NA,
              mcc_bootstrap_reps = bootstrap_reps,
              positive_class = positive_class,
              confidence = confidence,
              ts = ts
          )
      )
  }

  # --- Calculate Observed MCC and Chi-Square ---
  observed_mcc <- calculate_mcc_from_matrix(conf_matrix)
  chi_test <- chisq.test(conf_matrix, correct = TRUE)

  # --- Bootstrap Test for p-value ---
  bootstrap_p_val <- CI <- NA
  if (bootstrap_reps > 0) {
    if (ts) {
      # Time series bootstrap using the tsboot function from the boot package
      data_matrix <- cbind(as.character(actual), as.character(predicted))
      boot_out <- boot::tsboot(data_matrix, statistic = calculate_mcc, ordered_levels = ordered_levels,
                               R = bootstrap_reps, l = l, sim = sim, ...)
      mcc_boot_dist <- boot_out$t
    } else {
      # Standard bootstrap for non-time series data
      mcc_boot_dist <- replicate(bootstrap_reps, {
        i <- sample(seq_along(predicted), replace = TRUE)
        boot_matrix <- table(Actual = factor(c(actual[i], ordered_levels, rev(ordered_levels)), levels = ordered_levels),
                             Predicted = factor(c(predicted[i], ordered_levels, ordered_levels), levels = ordered_levels)) - 1
        calculate_mcc_from_matrix(boot_matrix)
      })
    }

    # Calculate the two-tailed p-value
    # https://github.com/mthulin/boot.pval/blob/main/R/boot.pval.R
    pval_precision <- 1 / bootstrap_reps
    alpha_seq <- seq(pval_precision, 1 - pval_precision, pval_precision)
    bounds <- cbind(quantile(mcc_boot_dist, probs = alpha_seq/2),
                    quantile(mcc_boot_dist, probs = 1 - alpha_seq/2))
    bootstrap_p_val <- alpha_seq[which.min(0 >= bounds[,1] & 0 <= bounds[,2])]

    # Calculate confidence intervals
    CI <- quantile(mcc_boot_dist, probs = c((1 - confidence) / 2, 1 - (1 - confidence) / 2))
  }

  # --- Output results ---
  list(
    confusion_matrix = conf_matrix,
    chi_square_test = chi_test,
    mcc = observed_mcc,
    mcc_bootstrap_pv = bootstrap_p_val,
    mcc_bootstrap_ci = CI,
    mcc_bootstrap_reps = bootstrap_reps,
    positive_class = positive_class,
    confidence = confidence,
    ts = ts
  )
}
