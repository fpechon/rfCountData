#' tuneRF
#' @name tuneRF
#' @description Starting with the default value of mtry, search for the optimal value 
#' (with respect to Out-of-Bag error estimate) of mtry for randomForest.
#' @usage tuneRF(x, offset, y, mtryStart, ntreeTry=50, stepFactor=2, improve=0.05,
#' trace=TRUE, plot=TRUE, doBest=FALSE, ...)
#' @param x, matrix or data frame of predictor variables,
#' @param offset, vector with log exposure,
#' @param y, response vector,
#' @param mtryStart, starting value of mtry; default is the same as in \code{rfPoisson}
#' @param ntreeTry, number of trees used at the tuning step
#' @param stepFactor, at each iteration, mtry is inflated (or deflated) by this value
#' @param improve, the (relative) improvement in OOB error must be by this  much for the search to continue
#' @param trace, whether to print the progress of the search
#' @param plot, whether to plot the OOB error as function of mtry
#' @param doBest, whether to run a forest using the optimal mtry found
#' @param ..., options to be given to \code{rfPoisson}
#' @return If \code{doBest=FALSE} (default), it returns a matrix whose first 
#' column contains the mtry values searched, and the second column the corresponding OOB error.
#' If \code{doBest=TRUE}, it returns the \code{rfPoisson}
#' object produced with the optimal \code{mtry}.
#' @export

tuneRF <- function(x, offset, y, mtryStart = floor(ncol(x)/3), ntreeTry=50, stepFactor=2,
                   improve=0.05, trace=TRUE, plot=TRUE, doBest=FALSE, ...) {
  if (improve < 0) stop ("improve must be non-negative.")
  errorOld <- rfPoisson(x, offset, y, mtry=mtryStart, ntree=ntreeTry,
                 keep.forest=FALSE, ...)$mse[ntreeTry]
  if (errorOld < 0) stop("Initial setting gave 0 error and no room for improvement.")
  if (trace) {
    cat("mtry =", mtryStart, " OOB error =",errorOld, "\n")
  }

  oobError <- list()
  oobError[[1]] <- errorOld
  names(oobError)[1] <- mtryStart  
  
  for (direction in c("left", "right")) {
    if (trace) cat("Searching", direction, "...\n")
    Improve <- 1.1*improve
    mtryBest <- mtryStart
    mtryCur <- mtryStart
    while (Improve >= improve) {
      mtryOld <- mtryCur
      mtryCur <- if (direction == "left") {
        max(1, ceiling(mtryCur / stepFactor))
      } else {
        min(ncol(x), floor(mtryCur * stepFactor))
      }
      if (mtryCur == mtryOld) break
      errorCur <- rfPoisson(x, offset, y, mtry=mtryCur, ntree=ntreeTry,
                     keep.forest=FALSE, ...)$mse[ntreeTry]
      if (trace) {
        cat("mtry =",mtryCur, "\tOOB error =", errorCur, "\n")
      }
      oobError[[as.character(mtryCur)]] <- errorCur
      Improve <- 1 - errorCur/errorOld
      cat(Improve, improve, "\n")
      if (Improve > improve) {
        errorOld <- errorCur
        mtryBest <- mtryCur
      }
    }
  }
  mtry <- sort(as.numeric(names(oobError)))
  res <- unlist(oobError[as.character(mtry)])
  res <- cbind(mtry=mtry, OOBError=res)

  if (plot) {
    plot(res, xlab=expression(m[try]), ylab="OOB Error", type="o", log="x",
         xaxt="n")
    axis(1, at=res[,"mtry"])
  }

  if (doBest) 
    res <- rfPoisson(x, offset, y, mtry=res[which.min(res[,2]), 1], ...)
  
  res
}
