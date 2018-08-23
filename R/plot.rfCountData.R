#' plot.rfCountData
#' @name plot.rfCountData
#' @param x, a rfCountData object
#' @param training, Should the Poisson deviance of the training set be plotted ? Defaults to \code{TRUE}.
#' @param testing, Should the Poisson deviance of the training set be plotted ? Defaults to \code{TRUE}, if a test set was used previously.
#' @export
plot.rfCountData <- function(x,  main=deparse(substitute(x)), training=TRUE, testing=TRUE,...) {
  if(x$type == "unsupervised")
    stop("No plot for unsupervised randomForest.")
  test <- !is.null(x$test$mse) & testing 
  d = data.frame(ntree = 1:x$ntree, training = x$mse, testing=NA)
  if (test)
    d$testing = x$test$mse
  
  p = ggplot2::ggplot(d)+ggplot2::xlab("trees") + ggplot2::ylab("Poisson Deviance")
  if (training)
  p = p + ggplot2::geom_line(ggplot2::aes(x = ntree, y = training, colour = "Training set")) 
  if(test) {
    d$testing = x$test$mse
    p = p + ggplot2::geom_line(ggplot2::aes(x = ntree, y = testing, colour="Testing set"))
  }
    p +ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "bottom")
}

  
