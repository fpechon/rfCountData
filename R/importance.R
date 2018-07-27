#' Extract variable importance measure
#' @param x, an object of class \code{\link{rfCountData}}.
#' @return A vector of importance measure, one item for each predictor variable.
#' @details The measure is computed from permuting OOB data: For each tree, the prediction error 
#' on the out-of-bag portion of the data is recorded (deviance). 
#' Then the same is done after permuting each predictor variable. 
#' The difference between the two are then averaged over all trees. 
#' @seealso \link{rfPoisson}, \link{varImpPlot}.
#' @export
importance <- function(x) {
    if (!inherits(x, "rfCountData"))
        stop("x is not of class rfCountData")
    hasImp <- !is.null(dim(x$importance)) || ncol(x$importance) == 2
    if (!hasImp)
        stop("That measure has not been computed")
    imp <- x$importance
    imp <- imp[, 1, drop=FALSE]
    imp
}
