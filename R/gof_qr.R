# https://github.com/vlyubchich/MML/blob/master/R/gof_qr.R

gof_qr <- function(obs, pred, quantiles = NULL) {
    pred <- cbind(pred)
    if (is.null(quantiles)) {
        # if not specified, try to get quantiles from the pred names
        quantiles <- colnames(pred)
        quantiles <- sapply(base::strsplit(quantiles, "[= ]+"), function(x) x[2])
        quantiles <- as.numeric(quantiles)
    }
    if (!is.null(quantiles) && (length(quantiles) != ncol(pred))) {
        stop("number of columns in 'pred' should match length(quantiles). Try setting 'quantiles = NULL'")
    }
    # pv function as in Eq 2 of Haupt et al. (2011)
    pv <- function(v, u) {
        (v - (u < 0)) * u
    }
    out <- sapply(1:ncol(pred), function(i) {
        e <- obs - pred[,i]
        # R1 function as in Eq 10 of Haupt et al. (2011)
        R1 <- 1 - sum(pv(quantiles[i], e)) / sum(pv(quantiles[i], obs - stats::quantile(obs, probs = quantiles[i])))
        # ATME function as in Eq 11 of Haupt et al. (2011)
        ATWE <- mean(pv(quantiles[i], e))
        c(R1 = R1, ATWE = ATWE)
    })
    colnames(out) <- paste0("quantile= ", quantiles)
    out
}
