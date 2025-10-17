#' Detect Rare Events Using Isolation Forest
#'
#' This function applies an isolation forest algorithm to detect rare events (anomalies)
#' in a dataset using the `isotree` package. It fits a model, computes anomaly scores,
#' and classifies observations based on a specified threshold.
#'
#' @param x A numeric matrix or data frame with no missing values (unless
#'   `missing_action` is set to handle them). Rows are observations, and columns are features.
#' @param ndim Integer, number of dimensions to split at each node (default: 1).
#' @param ntrees Integer, number of trees in the forest (default: 500).
#' @param missing_action Character, how to handle missing values: 'fail', 'impute', or 'auto' (default: 'fail').
#' @param scoring_metric Character, scoring method for anomalies: 'depth', 'adj_depth', or 'density' (default: 'adj_depth').
#' @param penalize_range Logical, whether to penalize features with smaller ranges (default: TRUE).
#' @param nthreads Integer, number of threads for parallel processing (default: parallel::detectCores() - 1).
#' @param threshold Numeric, anomaly score threshold for classifying rare events (default: 0.5).
#' @param sample_size Integer or fraction, number or proportion of rows to sample for training (default: NULL, uses all rows).
#' @param ... Additional arguments passed to `isotree::isolation.forest`.
#' @return A list containing:
#'   \itemize{
#'     \item \code{scores}: Numeric vector of anomaly scores for each observation.
#'     \item \code{is_anomaly}: Logical vector indicating whether each observation is an anomaly (score > threshold).
#'     \item \code{model}: The fitted `isotree` isolation forest model.
#'   }
#' @examples
#' \dontrun{
#' set.seed(123)
#' data <- matrix(rnorm(1000), nrow = 500)
#' data[1:5, ] <- data[1:5, ] + 10  # Add some outliers
#' result <- rare_iforest(data, ntrees = 100, threshold = 0.7)
#' table(result$is_anomaly)
#' plot(data, col = ifelse(result$is_anomaly, "red", "blue"), pch = 19)
#'
#' result <- rare_iforest(iris[, 1:4], threshold = 0.5)
#' table(result$is_anomaly)
#'
#' result <- rare_iforest(swiss, threshold = 0.5)
#' table(result$is_anomaly)
#' }
#' @importFrom isotree isolation.forest predict.isolation_forest
#' @export
#'
rare_iforest <- function(x,
                         ndim = 1,
                         ntrees = 500,
                         missing_action = "fail",
                         scoring_metric = "adj_depth",
                         penalize_range = TRUE,
                         nthreads = parallel::detectCores() - 1,
                         threshold = 0.5,
                         sample_size = NULL,
                         ...) {
    # Input validation
    if (!is.matrix(x) && !is.data.frame(x)) {
        stop("`x` must be a matrix or data frame.")
    }
    # if (is.data.frame(x) && !all(sapply(x, is.numeric))) {
    #     stop("All columns in `x` must be numeric.")
    # }
    if (missing_action == "fail" && any(is.na(x))) {
        stop("`x` contains missing values, and `missing_action` is set to 'fail'.")
    }
    if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
        stop("`threshold` must be a numeric value between 0 and 1.")
    }
    if (!requireNamespace("isotree", quietly = TRUE)) {
        stop("Package 'isotree' is required. Please install it using install.packages('isotree').")
    }

    # Handle sample_size for subsampling
    if (!is.null(sample_size)) {
        if (sample_size <= 1) {
            sample_size <- round(sample_size * nrow(x))
        }
        if (sample_size < 1 || sample_size > nrow(x)) {
            stop("`sample_size` must be a positive integer or fraction between 0 and 1.")
        }
        sample_idx <- sample(seq_len(nrow(x)), size = sample_size)
        x_sample <- x[sample_idx, , drop = FALSE]
    } else {
        x_sample <- x
    }

    # Fit isolation forest model
    model <- isotree::isolation.forest(
        x_sample,
        ndim = ndim,
        ntrees = ntrees,
        missing_action = missing_action,
        scoring_metric = scoring_metric,
        penalize_range = penalize_range,
        nthreads = max(1, as.integer(nthreads)),
        ...
    )

    # Predict anomaly scores
    scores <- predict(model, x)

    # Classify anomalies based on threshold
    is_anomaly <- scores > threshold

    # Return results
    return(list(
        scores = scores,
        is_anomaly = is_anomaly,
        model = model
    ))
}
