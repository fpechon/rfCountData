#' rfcv
#' @description This function shows the cross-validated prediction performance of models with sequentially reduced number of predictors (ranked by variable importance) via a nested cross-validation procedure.
#' @usage rfcv(trainx, trainy, cv.fold=5, scale="log", step=0.5, mtry=function(p) max(1, floor(sqrt(p))), recursive=FALSE, ...)
#' @param trainx, matrix or data frame containing columns of predictor variables
#' @param trainoffset, vector of offset, must have length equal to the number of rows in \code{trainx}
#' @param trainy, vector of response, must have length equal to the number of rows in \code{trainx}
#' @param cv.fold, number of folds in the cross-validation
#' @param scale, if \code{"log"}, reduce a fixed proportion (\code{step}) of variables at each step, otherwise reduce step variables at a time
#' @param step,	if \code{log=TRUE}, the fraction of variables to remove at each \code{step}, else remove this many variables at a time
#' @param mtry, a function of number of remaining predictor variables to use as the \code{mtry} parameter in the \code{rfPoisson} call
#' @param recursive, whether variable importance is (re-)assessed at each step of variable reduction
#' @param ..., other arguments passed on to \code{rfPoisson}
#' @author Andy Liaw
#' @references Svetnik, V., Liaw, A., Tong, C. and Wang, T., “Application of Breiman's Random Forest to Modeling Structure-Activity Relationships of Pharmaceutical Molecules”, MCS 2004, Roli, F. and Windeatt, T. (Eds.) pp. 334-343.
#' @export
rfcv <- function(trainx, trainoffset, trainy, cv.fold=5, scale="log", step=0.5,
                 mtry=function(p) max(1, floor(sqrt(p))), recursive=FALSE,
                 ...) {
    n <- nrow(trainx)
    p <- ncol(trainx)
    if (scale == "log") {
        k <- floor(log(p, base=1/step))
        n.var <- round(p * step^(0:(k-1)))
        same <- diff(n.var) == 0
        if (any(same)) n.var <- n.var[-which(same)]
        if (! 1 %in% n.var) n.var <- c(n.var, 1)
    } else {
        n.var <- seq(from=p, to=1, by=step)
    }
    k <- length(n.var)
    cv.pred <- vector(k, mode="list")
    for (i in 1:k) cv.pred[[i]] <- trainy
    ## Generate the indices of the splits
    ## For regression, bin the response into 5 bins and stratify.
		f <- factor(rep(1:5, length=length(trainy))[order(order(trainy))])
    nlvl <- table(f)
    idx <- numeric(n)
    for (i in 1:length(nlvl)) {
        idx[which(f == levels(f)[i])] <-  sample(rep(1:cv.fold, length=nlvl[i]))
    }

    for (i in 1:cv.fold) {
        ## cat(".")
        all.rf <- rfPoisson(trainx[idx != i, , drop=FALSE],
                               trainoffset[idx != i],
                               trainy[idx != i],
                               trainx[idx == i, , drop=FALSE],
                               trainy[idx == i],
                               trainoffset[idx == i],
                               mtry=mtry(p), importance=TRUE, ...)
        cv.pred[[1]][idx == i] <- all.rf$test$predicted
        impvar <- (1:p)[order(all.rf$importance[,1], decreasing=TRUE)]
        for (j in 2:k) {
            imp.idx <- impvar[1:n.var[j]]
            sub.rf <-
                rfPoisson(trainx[idx != i, imp.idx, drop=FALSE],
                          trainoffset[idx != i],
                          trainy[idx != i],
                          trainx[idx == i, imp.idx, drop=FALSE],
                          trainy[idx == i],
                          trainoffset[idx == i],
                          mtry=mtry(n.var[j]), importance=recursive, ...)
            cv.pred[[j]][idx == i] <- sub.rf$test$predicted
            ## For recursive selection, use importance measures from the sub-model.
            if (recursive) {
                impvar <-
                    (1:length(imp.idx))[order(sub.rf$importance[,1], decreasing=TRUE)]
      }
      NULL
    }
    NULL
  }
  ## cat("\n")
  error.cv <- sapply(cv.pred, function(x) mean(x - trainy*log(x+(x==0))))
  names(error.cv) <- names(cv.pred) <- n.var
  list(n.var=n.var, error.cv=error.cv, predicted=cv.pred)
}
