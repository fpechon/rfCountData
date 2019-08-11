#' Partial dependence plot
#' @name partialPlot
#' @description Partial dependence plot gives a graphical
#'  depiction of the marginal effect of a variable on the response variable.
#' @usage partialPlot(x, pred.data, x.var, offset,
#' w, plot = TRUE, n.pt = min(length(unique(pred.data[, xname])), 51),
#' rug = TRUE, xlab=deparse(substitute(x.var)), ylab="",
#' main=paste("Partial Dependence on", deparse(substitute(x.var))),...)
#' @param x, an object of class \code{rfCountData}, which contains a \code{forest} component.
#' @param pred.data, a data frame used for contructing the plot, usually the training data used to contruct the random forest.
#' @param offset, a vector, the corresponding log-exposures of pred.data.
#' @param x.var, name of the variable for which partial dependence is to be examined.
#' @param w, weights to be used in averaging; if not supplied, mean is not weighted
#' @param plot, whether the plot should be shown on the graphic device.
#' @param n.pt, if \code{x.var} is continuous, the number of points on the grid for evaluating partial dependence.
#' @param rug, whether to draw hash marks at the bottom of the plot indicating the deciles of \code{x.var}.
#' @param xlab, label for the x-axis.
#' @param ylab, label for the y-axis.
#' @param main, main title for the plot.
#' @param ..., other graphical parameters to be passed on to \code{plot} or \code{lines}.
#' @details The function being plotted is defined as: $\tilde{f}(x) = 1/n* sum_{i=1}^n f(x, x_{iC})$,
#' where x is the variable for which partial dependence is sought, and x_{iC} is the other variables in the data.
#' @return A list with two components: \code{x} and \code{y}, which are the values used in the plot. The \code{rfCountData} object must contain the forest component; i.e., created with \code{rfPoisson(..., keep.forest=TRUE)}. This function runs quite slow for large data sets.
#' @author Andy Liaw \email{andy_liaw@merck.com}
#' @references Friedman, J. (2001). Greedy function approximation: the gradient boosting machine, Ann. of Stat.
#' @seealso \link{rfPoisson}
#' @export
partialPlot <-
    function (x, pred.data, x.var, offset, w, plot=TRUE, n.pt = min(length(unique(pred.data[, xname])), 51),
              rug = TRUE, xlab=deparse(substitute(x.var)), ylab="",
              main=paste("Partial Dependence on", deparse(substitute(x.var))),
              ...)
{
    if (is.null(x$forest))
        stop("The randomForest object must contain the forest.\n")
    x.var <- substitute(x.var)
    xname <- if (is.character(x.var)) x.var else {
        if (is.name(x.var)) deparse(x.var) else {
            eval(x.var)
        }
    }
    xv <- pred.data[, xname]
    n <- nrow(pred.data)
    if (missing(w)) w <- rep(1, n)
    if (is.factor(xv) && !is.ordered(xv)) {
        x.pt <- levels(xv)
        y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
            x.data <- pred.data
            x.data[, xname] <- factor(rep(x.pt[i], n), levels = x.pt)
            y.pt[i] <- weighted.mean(predict(x, x.data,offset), w, na.rm=TRUE)
        }
        if (plot) p=ggplot2::ggplot() + ggplot2::geom_bar(ggplot2::aes(x = x.pt,y =y.pt, fill=x.pt), stat="identity")+
            ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(main) + ggplot2::theme(legend.position = "none")
    } else {
      
        if (is.ordered(xv)) xv <- as.numeric(xv)
        x.pt <- seq(min(xv), max(xv), length = n.pt)
        y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
            x.data <- pred.data
            x.data[, xname] <- rep(x.pt[i], n)
            y.pt[i] <- weighted.mean(predict(x, x.data, offset), w, na.rm=TRUE)
        }
        if (plot) p=ggplot2::ggplot() +  ggplot2::geom_line(ggplot2::aes(x=x.pt, y=y.pt)) +  
          ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +  ggplot2::ggtitle(main)
        if (rug && plot) {
            if (n.pt > 10) {
                p = p+ggplot2::geom_rug(ggplot2::aes(x=quantile(xv, seq(0.1, 0.9, by = 0.1))))
            } else {
              p = p+ggplot2::geom_rug(ggplot2::aes(x=unique(xv)))
            }
        }
    }
    p
}
