#' @export
"print.rfCountData" <-
function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n")
  cat("               Type of random forest: ", x$type, "\n", sep="")
  cat("                     Number of trees: ", x$ntree, "\n",sep="")
  cat("No. of variables tried at each split: ", x$mtry, "\n\n", sep="")
  if(!is.null(x$mse)) {
    cat("           Poisson Loss Function: ", x$mse[length(x$mse)],
        "\n", sep="")
    if(!is.null(x$test$mse)) {
      cat("Test set Poisson Loss Function: ",
          round(x$test$mse[length(x$test$mse)], digits=2), "\n", sep="")
    }      
  }
}
