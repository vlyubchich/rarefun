#' Calculate and Test the Matthews Correlation Coefficient (MCC)
#'
#' Based on two vectors of binary values, compute the confusion matrix and the MCC
#' (mean square contingency coefficient).
#' The function implements two methods for testing significance:
#' 1. A parametric chi-square test.
#' 2. A non-parametric bootstrap test for a more robust p-value.
#'
#' @param x A vector of binary labels. Can be numeric (0/1),
#'   logical (TRUE/FALSE), character, or factor.
#' @param y A vector of binary labels, of the same type and
#'   length as \code{x}.
#' @param positive_class An optional value explicitly specifying the "positive"
#'   class label. If \code{NULL}, the function will infer it.
#' @param bootstrap_reps The number of bootstrap replicates for p-value calculation.
#'   Default is 999. Set to 0 to disable bootstrapping.
#' @param ts A logical flag indicating if the data are time series. Default is FALSE. If TRUE,
#'  a version of bootstrap for time series is applied to account for potential autocorrelation;
#'  the data are then assumed to be ordered in time.
#' @param confidence The confidence level for the bootstrap confidence interval (default: 0.95).
#'  Must be between 0 and 1 (exclusive).
#' @param l The block length for the time series bootstrap (default: 5). Only used if \code{ts} is TRUE,
#'  see \code{?boot::tsboot}.
#' @param sim The type of simulation for the time series bootstrap (default: "fixed").
#'  Options are "fixed" (moving block bootstrap) or "geom" (stationary bootstrap), see \code{?boot::tsboot}.
#' @param ... Additional arguments passed to \code{boot::tsboot} when \code{ts} is TRUE.
#'
#'
#' @return A list containing:
#'   \itemize{
#'   \item \code{confusion_matrix}: The 2x2 confusion matrix.
#'   \item \code{mcc}: The Matthews Correlation Coefficient (-1 to +1).
#'   \item \code{chi_square_test}: The output of \code{chisq.test()}. \code{$p.value} is the parametric p-value.
#'   \item \code{mcc_bootstrap_pv}: The p-value from the bootstrap permutation test.
#'   \item \code{mcc_bootstrap_ci}: The bootstrap confidence interval for the MCC.
#'   \item \code{mcc_bootstrap_reps}: The number of bootstrap replicates used.
#'   \item \code{positive_class}: The positive class label used.
#'   \item \code{confidence}: The confidence level for the bootstrap confidence interval.
#'   \item \code{ts}: Whether the data were treated as time series.
#'   }
#'
#' @details
#' The MCC is a robust metric for binary classification, especially on imbalanced data.
#' It ranges from -1 (perfect negative correlation) to +1 (perfect positive correlation).
#' For 2x2 confusion matrices:
#' \deqn{MCC = \frac{TP\,TN - FP\,FN}{\sqrt{(TP + FP)(TP + FN)(TN + FP)(TN + FN)}}}{MCC = (TP*TN - FP*FN)/sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))}
#' where:
#' - TP: True Positives
#' - TN: True Negatives
#' - FP: False Positives
#' - FN: False Negatives
#'
#' and
#' \deqn{|MCC| = \sqrt{\frac{\chi^2}{n}}}{|MCC| = sqrt(chi^2 / n)}
#'
#' The **chi-square test** is fast but assumes that each observation is independent,
#' which is often not true for time series data.
#'
#' **For Time Series:** Both standard chi-square and standard bootstrapping assume
#' data independence. If the \code{ts} parameter is set to TRUE, the function applies a
#' version of bootstrapping suitable for time series data to account for potential
#' autocorrelation. The data are assumed to be ordered in time. The block length \code{l}
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
#' The function uses the \code{boot} package for bootstrapping.
#'
#' @seealso \code{\link[stats]{chisq.test}}, \code{\link[boot]{tsboot}}
#'
#' @examples
#' # Example 1: A clear, significant correlation
#' x_vals <- rep(c(1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1), 3)
#' y_vals <- rep(c(1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0), 3)
#' mcc_results <- mcc(x_vals, y_vals, positive_class = 1)
#' print(mcc_results$confusion_matrix)
#' print(paste("MCC:", round(mcc_results$mcc, 3)))
#' print(paste("Chi-Square p-value:", round(mcc_results$chi_square_test$p.value, 4)))
#' print(paste("Bootstrap p-value:", round(mcc_results$mcc_bootstrap_pv, 4)))
#'
#' # Example 2: No significant correlation
#' x_rand <- rep(c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0), 3)
#' y_rand <- rep(c(1, 1, 0, 0, 1, 1, 0, 0, 1, 0), 3)
#' mcc_rand <- mcc(x_rand, y_rand, bootstrap_reps = 999)
#' print(paste("MCC:", round(mcc_rand$mcc, 3)))
#' print(paste("Bootstrap p-value:", round(mcc_rand$mcc_bootstrap_pv, 4)))
#'
#' # Example 3: Time series data with autocorrelation in both series
#' n <- 100
#' set.seed(12345)
#' ts_x <- rbinom(n, 1, 0.3)
#' ts_y <- ifelse(runif(n) < 0.3,
#'                ts_x,
#'                1 - ts_x)
#' # Introduce autocorrelation
#' for (i in 2:n) {
#'    if (runif(1) < 0.7) {
#'        ts_x[i] <- ts_x[i - 1]
#'        ts_y[i] <- ts_y[i - 1]
#'    }
#' }
#'
#' # Check autocorrelation at lag 1
#' par(mfrow = c(1, 2))
#' mcc(ts_x, dplyr::lag(ts_x), ts = TRUE, l = 7)$mcc
#' acf(ts_x, main = "ACF of X Time Series, treated as numeric")
#' mcc(ts_y, dplyr::lag(ts_y), ts = TRUE, l = 7)$mcc
#' acf(ts_y, main = "ACF of Y Time Series")
#'
#' # Visualize the time series
#' plot.ts(cbind(ts_x, ts_y),
#'        main = "Time Series of X and Y", xlab = "Time")
#'
#' # Calculate MCC treating data as time series
#' mcc_ts <- mcc(ts_x,
#'               ts_y,
#'               positive_class = 1,
#'               ts = TRUE, l = 7, bootstrap_reps = 999)
#' print(paste("MCC:", round(mcc_ts$mcc, 3)))
#' print(paste("Chi-Square p-value:", round(mcc_ts$chi_square_test$p.value, 4)))
#' print(paste("Bootstrap p-value:", round(mcc_ts$mcc_bootstrap_pv, 4)))
#' print(paste("Bootstrap CI:",
#'    paste(round(mcc_ts$mcc_bootstrap_ci, 3), collapse = " to ")))
#'
#' @importFrom stats chisq.test quantile rmultinom
#' @importFrom boot tsboot
#' @export
#'
mcc <- function(x,
                y,
                positive_class = NULL,
                bootstrap_reps = 999,
                confidence = 0.95,
                ts = FALSE,
                l = 5,
                sim = "fixed", ...) {

  # --- Capture argument names for confusion matrix labels ---
  x_name <- deparse(substitute(x))
  y_name <- deparse(substitute(y))

  # --- Internal helper function to compute MCC from counts ---
  calculate_mcc_from_counts <- function(counts) {
    TP <- as.numeric(counts[1])
    FN <- as.numeric(counts[2])
    FP <- as.numeric(counts[3])
    TN <- as.numeric(counts[4])
    numerator <- TP * TN - FP * FN
    denominator_parts <- c(TP + FP, TP + FN, TN + FP, TN + FN)
    if (any(denominator_parts == 0)) {
      # If any part of the denominator is zero, return 0 to avoid division by zero
      return(0)
    }
    denominator <- sqrt(prod(denominator_parts))
    numerator / denominator
  }

  # --- Internal helper functions for compact binary encoding ---
  encode_pair_codes <- function(x_binary, y_binary) {
    4L - (2L * as.integer(x_binary)) - as.integer(y_binary)
  }

  build_confusion_matrix <- function(counts, ordered_levels, xname, yname) {
    dimnames_list <- list(ordered_levels, ordered_levels)
    names(dimnames_list) <- c(xname, yname)
    matrix(counts, nrow = 2, byrow = TRUE, dimnames = dimnames_list)
  }

  calculate_bootstrap_mcc <- function(data) {
    counts <- tabulate(as.integer(data), nbins = 4L)
    calculate_mcc_from_counts(counts)
  }

  # --- Input Validation ---
  if (length(x) != length(y)) {
    stop("Input vectors 'x' and 'y' must have the same length.")
  }
  if (length(x) == 0) {
    stop("Input vectors must not be empty.")
  }
  unique_values <- sort(unique(c(as.character(x), as.character(y))))
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
  if (length(unique(x)) > 2) {
    stop("'x' vector must be binary (contain at most two unique values).")
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
    if (is.logical(x)) {
      positive_class <- TRUE
    } else if (is.numeric(x) && all(x %in% c(0, 1))) {
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
  x_binary <- as.character(x) == positive_class
  y_binary <- as.character(y) == positive_class
  pair_codes <- encode_pair_codes(x_binary, y_binary)
  conf_counts <- tabulate(pair_codes, nbins = 4L)
  conf_matrix <- build_confusion_matrix(conf_counts, ordered_levels, x_name, y_name)

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
  observed_mcc <- calculate_mcc_from_counts(conf_counts)
  chi_test <- chisq.test(conf_matrix, correct = TRUE)

  # --- Bootstrap Test for p-value ---
  bootstrap_p_val <- CI <- NA
  if (bootstrap_reps > 0) {
    if (ts) {
      # Time series bootstrap using the tsboot function from the boot package
      boot_out <- boot::tsboot(pair_codes, statistic = calculate_bootstrap_mcc,
                               R = bootstrap_reps, l = l, sim = sim, ...)
      mcc_boot_dist <- boot_out$t
    } else {
      # Standard bootstrap for non-time series data: resample the 4 pair categories directly.
      boot_counts <- rmultinom(bootstrap_reps,
                               size = length(pair_codes),
                               prob = conf_counts / sum(conf_counts))
      TP <- boot_counts[1, ]
      FN <- boot_counts[2, ]
      FP <- boot_counts[3, ]
      TN <- boot_counts[4, ]
      denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
      mcc_boot_dist <- ifelse(denominator == 0,
                              0,
                              (TP * TN - FP * FN) / denominator)
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
