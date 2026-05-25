#' Partial Dependence Plot for Quantile Random Forest
#'
#' Get a data frame for partial dependence plots from a quantile random forest.
#' The plots can be obtained using plotting functions from other packages.
#'
#' @param object a \code{\link[ranger]{ranger}} quantile random forest object.
#' @param pred.var character string giving the names of the predictor variables of
#' interest (see \code{\link[pdp]{partial}}).
#' @param Q a numeric vector of probabilities for which the plot is desired
#'   (default: \code{c(0.05, 0.5, 0.95))}. Values should be between 0 and 1.
#' @param ... other arguments passed to \code{\link[pdp]{partial}}.
#'
#' @return A \code{data.frame} with the following columns:
#' the predictor variables supplied in \code{pred.var},
#' response variable (the name is extracted from the \code{ranger} function call),
#' and, if \code{length(Q) > 1}, one more column named \code{"Quantile"}
#' (for an appropriate order of labels in the plots, \code{Quantile} is a factor).
#'
#' @seealso \code{\link[pdp]{partial}}, \code{\link[ranger]{ranger}}
#'
#' @keywords quantile forest partial
#'
#' @importFrom pdp partial
#' @importFrom stats predict
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ranger", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'
#'   # fit a quantile random forest
#'   qrf <- ranger::ranger(Examination ~ ., data = swiss, quantreg = TRUE, num.trees = 10)
#'
#'   # 1) a plot for one quantile
#'   partial_qrf(qrf, pred.var = "Agriculture", Q = 0.5) |>
#'       ggplot(aes(Agriculture, Examination)) +
#'           geom_line()
#'
#'   # 2) a plot for several quantiles, with color scale
#'   if (requireNamespace("viridis", quietly = TRUE)) {
#'     library(viridis)
#'     partial_qrf(qrf, pred.var = "Agriculture") |>
#'        ggplot(aes(Agriculture, Examination, group = Quantile, color = Quantile)) +
#'            geom_line() +
#'            scale_color_viridis(discrete = TRUE) +
#'            theme_minimal()
#'   }
#'
#'   # 3) a plot for one quantile for 2 predictors (use small grid.resolution to save time)
#'   df <- partial_qrf(qrf, pred.var = c("Agriculture", "Catholic"), Q = 0.5,
#'       grid.resolution = 5)
#'   if (requireNamespace("viridis", quietly = TRUE)) {
#'     ggplot(df, aes(Agriculture, Catholic)) +
#'         geom_tile(aes(fill = Examination)) +
#'         scale_fill_viridis() +
#'         theme_minimal() +
#'         labs(fill = "Median\nexamination")
#'   }
#'
#'   # 4) a plot for each predictor
#'   varnames <- qrf$forest$independent.variable.names
#'   qs <- c(0.05, 0.5, 0.95)
#'   ddf <- lapply(varnames, function(vn) partial_qrf(qrf, pred.var = vn, Q = qs))
#'   if (requireNamespace("viridis", quietly = TRUE) &&
#'       requireNamespace("patchwork", quietly = TRUE)) {
#'     library(viridis)
#'     library(patchwork)
#'     yrange <- range(sapply(ddf, function(x) range(x[, 2])))
#'     plist <- lapply(seq_along(varnames), function(vi) {
#'         ggplot(ddf[[vi]], aes(x = .data[[varnames[vi]]], y = Examination,
#'             group = Quantile, color = Quantile)) +
#'             geom_line() +
#'             scale_color_viridis(discrete = TRUE) +
#'             theme_minimal() +
#'             ylim(yrange[1], yrange[2]) +
#'             theme(axis.title.y = element_blank())
#'     })
#'     wrap_plots(plist) + plot_layout(guides = "collect")
#'   }
#' }
#'}
#'
partial_qrf <- function(object, pred.var
                        ,Q = c(0.05, 0.5, 0.95)
                        ,...)
{
    Q <- sort(unique(Q))
    predfun <- function(object, newdata) {
        qpred <- predict(object, newdata, type = "quantiles", quantiles = Q)$predictions
        tmp <- c(apply(qpred, 2, mean))
        qq <- paste0("q", Q)
        # qq <<-
        qq <- gsub("\\.", "_", qq)
        names(tmp) <- qq
        return(tmp)
    }
    pdpout <- pdp::partial(object
                           ,pred.var = pred.var
                           ,pred.fun = predfun
                           ,...)
    if ("call" %in% ls(object)) {
        response <- as.character(object$call)
        response <- base::strsplit(response[2], "\\ ")[[1]][1]
        colnames(pdpout)[length(pred.var) + 1] <- response
    }
    if (length(Q) > 1) {
        colnames(pdpout)[ncol(pdpout)] <- "Quantile"
        # pdpout$Quantile <- factor(pdpout$Quantile, levels = rev(qq))
    }
    return(pdpout)
}
