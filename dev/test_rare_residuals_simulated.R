rm(list = ls())

# Simulated benchmarking script for rare_residuals()
# - Evaluates detection performance on non-seasonal and seasonal synthetic data
# - Compares Isolation Forest and DBSCAN outputs from rare_residuals(method = "all")
# - Reports precision/recall/F1/balanced accuracy/MCC and simple pass thresholds

if (!"package:rarefun" %in% search()) {
    if (requireNamespace("devtools", quietly = TRUE)) {
        devtools::load_all(".", quiet = TRUE)
    } else {
        stop("Load rarefun first (library(rarefun)) or install devtools to use load_all('.').")
    }
}

calc_metrics <- function(truth, pred) {
    truth <- as.logical(truth)
    pred <- as.logical(pred)

    tp <- sum(pred & truth)
    fp <- sum(pred & !truth)
    tn <- sum(!pred & !truth)
    fn <- sum(!pred & truth)

    safe_div <- function(num, den) if (den == 0) NA_real_ else num / den

    precision <- safe_div(tp, tp + fp)
    recall <- safe_div(tp, tp + fn)
    specificity <- safe_div(tn, tn + fp)
    accuracy <- safe_div(tp + tn, tp + tn + fp + fn)
    f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
        NA_real_
    } else {
        2 * precision * recall / (precision + recall)
    }
    bal_acc <- mean(c(recall, specificity), na.rm = TRUE)

    mcc_num <- (tp * tn) - (fp * fn)
    mcc_den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc <- if (mcc_den == 0) NA_real_ else mcc_num / mcc_den

    data.frame(
        tp = tp,
        fp = fp,
        tn = tn,
        fn = fn,
        precision = precision,
        recall = recall,
        specificity = specificity,
        f1 = f1,
        accuracy = accuracy,
        balanced_accuracy = bal_acc,
        mcc = mcc
    )
}

simulate_nonseasonal <- function(n = 900,
                                 event_rate = 0.02,
                                 noise_sd = 1,
                                 event_size = 6) {
    time <- seq_len(n)
    baseline <- 5 + 0.01 * time + 0.6 * sin(2 * pi * time / 120)
    value <- baseline + stats::rnorm(n, sd = noise_sd)

    n_events <- max(1, floor(n * event_rate))
    event_idx <- sort(sample(20:(n - 20), n_events, replace = FALSE))
    event_sign <- sample(c(-1, 1), size = n_events, replace = TRUE)
    value[event_idx] <- value[event_idx] + event_sign * event_size * noise_sd

    list(
        x = value,
        event_idx = event_idx,
        truth = seq_len(n) %in% event_idx
    )
}

simulate_seasonal <- function(n = 10 * 365,
                              period = 365.25,
                              event_rate = 0.01,
                              noise_sd = 1,
                              event_size = 6,
                              start_date = as.Date("2000-01-01")) {
    time <- seq_len(n)
    date <- seq(start_date, by = "day", length.out = n)

    baseline <-
        12 +
        0.001 * time +
        8 * sin(2 * pi * time / period) +
        2 * cos(4 * pi * time / period)

    value <- baseline + stats::rnorm(n, sd = noise_sd)

    n_events <- max(1, floor(n * event_rate))
    event_idx <- sort(sample(30:(n - 30), n_events, replace = FALSE))
    event_sign <- sample(c(-1, 1), size = n_events, replace = TRUE)
    value[event_idx] <- value[event_idx] + event_sign * event_size * noise_sd

    list(
        x = data.frame(time = date, value = value),
        event_idx = event_idx,
        truth = seq_len(n) %in% event_idx
    )
}

evaluate_iforest_tuned <- function(truth, score_iforest,
                                   thresholds = seq(0.35, 0.90, by = 0.05)) {
    scored <- lapply(thresholds, function(th) {
        pred <- score_iforest > th
        met <- calc_metrics(truth, pred)
        met$threshold <- th
        met
    })
    scored <- do.call(rbind, scored)
    scored <- scored[order(scored$f1, scored$balanced_accuracy, decreasing = TRUE), ]
    scored[1, , drop = FALSE]
}

run_one_rep <- function(seasonal = FALSE,
                        iforest_args = list(ntrees = 300, threshold = 0.60),
                        dbscan_args = list(minPts = 5)) {
    sim <- if (seasonal) {
        simulate_seasonal()
    } else {
        simulate_nonseasonal()
    }

    fit <- rare_residuals(
        sim$x,
        seasonal = seasonal,
        period = if (seasonal) 365.25 else NULL,
        method = "all",
        iforest_args = iforest_args,
        dbscan_args = dbscan_args
    )

    d <- fit$data
    truth <- sim$truth

    iforest_default <- calc_metrics(truth, d$is_anomaly_iforest)
    iforest_default$detector <- "iforest_default"

    iforest_tuned <- evaluate_iforest_tuned(truth, d$score_iforest)
    iforest_tuned$detector <- "iforest_tuned"

    dbscan <- calc_metrics(truth, d$is_anomaly_dbscan)
    dbscan$detector <- "dbscan"

    out <- rbind(iforest_default, iforest_tuned, dbscan)
    out
}

summarize_results <- function(res_df) {
    metrics <- c("precision", "recall", "f1", "balanced_accuracy", "mcc")
    split_df <- split(res_df, res_df$detector)

    by_detector <- lapply(names(split_df), function(det_name) {
        dat <- split_df[[det_name]]
        row <- data.frame(detector = det_name)
        for (m in metrics) {
            row[[paste0(m, "_mean")]] <- mean(dat[[m]], na.rm = TRUE)
            row[[paste0(m, "_sd")]] <- stats::sd(dat[[m]], na.rm = TRUE)
        }
        row
    })

    do.call(rbind, by_detector)
}

run_benchmark <- function(n_reps = 40, seasonal = FALSE, seed = 20260303) {
    set.seed(seed)
    out <- lapply(seq_len(n_reps), function(i) {
        rep_res <- run_one_rep(seasonal = seasonal)
        rep_res$replicate <- i
        rep_res$series_type <- if (seasonal) "seasonal" else "nonseasonal"
        rep_res
    })
    out <- do.call(rbind, out)
    list(per_rep = out, summary = summarize_results(out))
}

reps_nonseasonal <- as.integer(Sys.getenv("RAREFUN_SIM_REPS_NONSEASONAL", "40"))
reps_seasonal <- as.integer(Sys.getenv("RAREFUN_SIM_REPS_SEASONAL", "40"))

cat("\nRunning simulation benchmark for rare_residuals()...\n")
cat(sprintf("Replicates -> nonseasonal: %d, seasonal: %d\n", reps_nonseasonal, reps_seasonal))
nonseasonal_res <- run_benchmark(n_reps = reps_nonseasonal, seasonal = FALSE, seed = 20260303)
seasonal_res <- run_benchmark(n_reps = reps_seasonal, seasonal = TRUE, seed = 20260304)

summary_all <- rbind(nonseasonal_res$summary, seasonal_res$summary)

cat("\n=== Mean performance across replicates ===\n")
print(summary_all)

# Simple capability checks (adjust thresholds as needed for your use case)
check_capability <- function(df, min_recall = 0.50, min_bal_acc = 0.75) {
    recall_col <- "recall_mean"
    bal_col <- "balanced_accuracy_mean"

    df$capable <- with(df, get(recall_col) >= min_recall & get(bal_col) >= min_bal_acc)
    df
}

checks <- check_capability(summary_all)
cat("\n=== Capability checks (recall >= 0.50 and balanced accuracy >= 0.75) ===\n")
print(checks[, c("detector", "recall_mean", "balanced_accuracy_mean", "capable")])

# Save results for later inspection
dir.create("dev/results", showWarnings = FALSE, recursive = TRUE)
write.csv(summary_all,
          file = "dev/results/rare_residuals_simulated_summary.csv",
          row.names = FALSE)
write.csv(rbind(nonseasonal_res$per_rep, seasonal_res$per_rep),
          file = "dev/results/rare_residuals_simulated_per_rep.csv",
          row.names = FALSE)

cat("\nSaved:\n")
cat(" - dev/results/rare_residuals_simulated_summary.csv\n")
cat(" - dev/results/rare_residuals_simulated_per_rep.csv\n")
