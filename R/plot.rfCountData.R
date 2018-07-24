#' plot.rfCountData
#' @name plot.rfCountData
#' @export
plot.rfCountData <- function(x, type="l", main=deparse(substitute(x)), ...) {
  if(x$type == "unsupervised")
    stop("No plot for unsupervised randomForest.")
  test <- !is.null(x$test$mse) 
  err <- x$mse
  if(test) err <- cbind(err, x$test$mse)

  if(test) {
    colnames(err) <- c("OOB", "Test")
    matplot(1:x$ntree, err, type = type, xlab="trees", ylab="Error",
            main=main, ...)
  } else {
    matplot(1:x$ntree, err, type = type, xlab="trees", ylab="Error",
            main=main, ...)
  }
  invisible(err)
}

  
