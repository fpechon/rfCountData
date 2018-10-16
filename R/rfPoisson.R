## mylevels() returns levels if given a factor, otherwise 0.
mylevels <- function(x) if (is.factor(x)) levels(x) else 0


#' rfPoisson
#' @name rfPoisson
#' @description \code{rfPoisson} implements Breiman's random forest algorithm (based on 
#' Breiman and Cutler's original Fortran code) for regression and has been modified 
#' to be used with Poisson data that have different observation periods.\cr
#' More specifically, the best split is the one that will maximise the decrease of the poisson deviance. An offset
#' has also been introduced to accomodate for different times of exposure. The offset the log of the exposure.
#' @param y a vector of Poisson responses
#' @param x a data frame or a matrix of predictors.
#' @param offset a vector of same size as y, corresponding to the log of observation time (e.g. log of exposure). Default is 0.
#' @param xtest a data frame or matrix (like \code{x}) containing predictors for the test set.
#' @param ytest response for the test set
#' @param offsettest Offset for the test set (like \code{offset})
#' @param ntree, Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. 
#' @param mtry Number of variables randomly sampled as candidates at each split.  Default is p/3, where p is the number of variables in \code{x}.
#' @param replace Should sampling of cases be done with or without replacement?
#' @param sampsize Size(s) of sample to draw.
#' @param nodesize Minimum size of terminal nodes.  Setting this number 
#' larger causes smaller trees to be grown (and thus take less time). Default is 5000.
#' @param maxnodes Maximum number of terminal nodes trees in the forest can have.  
#' If not given, trees are grown to the maximum possible 
#' (subject to limits by \code{nodesize}).  If set larger than maximum possible, a warning is issued.
#' @param importance Should importance of predictors be assessed? Default is \code{TRUE}.
#' @param nPerm Not yet implemented. Number of times the OOB data are permuted per tree for 
#' assessing variable importance.  Number larger than 1 gives slightly more stable estimate, 
#' but not very effective.
#' @param do.trace If set to \code{TRUE}, give a more verbose output as \code{rfPoisson} is run.
#'   If set to some integer, then running output is printed for every \code{do.trace} trees.
#' @param keep.forest If set to \code{FALSE}, the forest will not be retained in the output object.
#'   If \code{xtest} is given, defaults to \code{FALSE}.
#' @param keep.inbag Should an \code{n} by \code{ntree} matrix be returned that keeps track of which 
#' samples are `in-bag' in which trees (but not how many times, if sampling with replacement
#' @param ... other parameters passed to lower functions.
#' @return TBC
#' @note TBC
#' @references 
#' R package \emph{randomForest},  \url{https://cran.r-project.org/package=randomForest}\cr
#' Breiman, L. (2001), \emph{Random Forests}, Machine Learning 45(1), 5-32.\cr
#' Breiman, L (2002), ``Manual On Setting Up, Using, And Understanding Random Forests V3.1'', \url{https://www.stat.berkeley.edu/~breiman/Using_random_forests_V3.1.pdf}.
#' @examples 
#' if (!require(CASdatasets)) install.packages("CASdatasets", repos = "http://cas.uqam.ca/pub/R/", type="source")
#' require(CASdatasets)
#' data("freMTPLfreq")
#' library(rfCountData)
#' m0 = rfPoisson(y = freMTPLfreq[1:10000,]$ClaimNb,
#'                   offset = log(freMTPLfreq[1:10000,]$Exposure),
#'                   x = freMTPLfreq[1:10000,c("Region", "Power", "DriverAge")],
#'                   ntree = 20)
#' predict(m0, newdata = freMTPLfreq[10001:10050,c("Region", "Power", "DriverAge")], 
#' offset = log(freMTPLfreq[10001:10050,"Exposure"]))
#' @author Florian Pechon, \email{florian.pechon@uclouvain.be}, based on the package randomForest by
#'  Andy Liaw \email{andy_liaw@merck.com} and Matthew Wiener \email{matthew_wiener@merck.com} based 
#'  on original Fortran code by Leo Breiman and Adele Cutler.
#' @export

rfPoisson <-
    function(x, offset=NULL, y=NULL,  xtest=NULL, ytest=NULL, offsettest=NULL, ntree=500,
             mtry=max(floor(ncol(x)/3), 1),
             replace=TRUE,
             sampsize = if (replace) nrow(x) else ceiling(.632*nrow(x)),
             nodesize = 5000,
             maxnodes=NULL,
             importance=TRUE, 
             nPerm=1,
             do.trace=FALSE,
             keep.forest=!is.null(y) && is.null(xtest),
             keep.inbag=FALSE, ...) {
    n <- nrow(x)
    p <- ncol(x)
    if (n == 0) stop("data (x) has 0 rows")
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)

    ## overcome R's lazy evaluation:
    keep.forest <- keep.forest

    testdat <- !is.null(xtest)
    if (testdat) {
        if (ncol(x) != ncol(xtest))
            stop("x and xtest must have same number of columns")
        ntest <- nrow(xtest)
        xts.row.names <- rownames(xtest)
    }

    ## Make sure mtry is in reasonable range.
    if (mtry < 1 || mtry > p)
        warning("invalid mtry: reset to within valid range")
    mtry <- max(1, min(p, round(mtry)))
        if (length(y) != n) stop("length of response must be the same as predictors")


    ## Check for NAs.
    if (any(is.na(x))) stop("NA not permitted in predictors")
    if (testdat && any(is.na(xtest))) stop("NA not permitted in xtest")
    if (any(is.na(y))) stop("NA not permitted in response")
    if (!is.null(ytest) && any(is.na(ytest))) stop("NA not permitted in ytest")

    if (is.data.frame(x)) {
        xlevels <- lapply(x, mylevels)
        ncat <- sapply(xlevels, length)
        ## Treat ordered factors as numerics.
        ncat <- ifelse(sapply(x, is.ordered), 1, ncat)
        x <- data.matrix(x)
        if(testdat) {
            if(!is.data.frame(xtest))
                stop("xtest must be data frame if x is")
            xfactor <- which(sapply(xtest, is.factor))
            if (length(xfactor) > 0) {
                for (i in xfactor) {
                    if (any(! levels(xtest[[i]]) %in% xlevels[[i]]))
                        stop("New factor levels in xtest not present in x")
                    xtest[[i]] <-
                        factor(xlevels[[i]][match(xtest[[i]], xlevels[[i]])],
                               levels=xlevels[[i]])
                }
            }
            xtest <- data.matrix(xtest)
        }
    } else {
        ncat <- rep(1, p)
		names(ncat) <- colnames(x)
        xlevels <- as.list(rep(0, p))
    }
    maxcat <- max(ncat)
    if (maxcat > 53)
        warning("Can not handle categorical predictors with more than 53 categories.")
    if (importance) {
        if (nPerm < 1) nPerm <- as.integer(1) else nPerm <- as.integer(nPerm)
            impout <- matrix(0.0, p, 2)
            impSD <- double(p)
            names(impSD) <- x.col.names
    } else {
        impout <- double(p)
        impSD <- double(1)
    }

    nsample <-  n
    Stratify <- length(sampsize) > 1
    if (Stratify) stop("sampsize should be of length one")
    ## For regression trees, need to do this to get maximal trees.
    nrnodes <- 2 * trunc(sampsize/max(1, nodesize - 4)) + 1
    if (!is.null(maxnodes)) {
        ## convert # of terminal nodes to total # of nodes
        maxnodes <- 2 * maxnodes - 1
        if (maxnodes > nrnodes) warning("maxnodes exceeds its max value.")
        nrnodes <- min(c(nrnodes, max(c(maxnodes, 1))))
    }
    ## Compiled code expects variables in rows and observations in columns.
    x <- t(x)
    storage.mode(x) <- "double"
    if (testdat) {
        xtest <- t(xtest)
        storage.mode(xtest) <- "double"
        if (is.null(ytest)) {
            ytest <- labelts <- 0
        } else {
            labelts <- TRUE
        }
    } else {
        xtest <- double(1)
        ytest <- double(1)
        ntest <- 1
        labelts <- FALSE
    }
    nt <- if (keep.forest) ntree else 1


      rfout <- .C("regRF",
                    data.matrix(x),
                    as.double(offset),
                    as.double(y),
                    as.integer(c(n, p)),
                    as.integer(sampsize),
                    as.integer(nodesize),
                    as.integer(nrnodes),
                    as.integer(ntree),
                    as.integer(mtry),
                    as.integer(c(importance, 
                                 nPerm)),
                    as.integer(ncat),
                    as.integer(maxcat),
                    as.integer(do.trace),
                    ypred = double(n),
                    impout = impout,
                    ndbigtree = integer(ntree),
                    nodestatus = matrix(integer(nrnodes * nt), ncol=nt),
                    leftDaughter = matrix(integer(nrnodes * nt), ncol=nt),
                    rightDaughter = matrix(integer(nrnodes * nt), ncol=nt),
                    nodepred = matrix(double(nrnodes * nt), ncol=nt),
                    bestvar = matrix(integer(nrnodes * nt), ncol=nt),
                    xbestsplit = matrix(double(nrnodes * nt), ncol=nt),
                    dev = double(ntree),
                    keep = as.integer(c(keep.forest, keep.inbag)),
                    replace = as.integer(replace),
                    testdat = as.integer(testdat),
                    xts = data.matrix(xtest),
                    ntest = as.integer(ntest),
                    yts = as.double(ytest),
                    offsetts = as.double(offsettest),
                    labelts = as.integer(labelts),
                    ytestpred = double(ntest),
                    devts = double(if (labelts) ntree else 1),
                    coef = double(2),
                    oob.times = integer(n),
                    inbag = if (keep.inbag)
                    matrix(integer(n * ntree), n) else integer(1),
                    DUP=FALSE,
                    PACKAGE="rfCountData")[c(14:24, 31:37)]
        ## Format the forest component, if present.
        if (keep.forest) {
            max.nodes <- max(rfout$ndbigtree)
            rfout$nodestatus <-
                rfout$nodestatus[1:max.nodes, , drop=FALSE]
            rfout$bestvar <-
                rfout$bestvar[1:max.nodes, , drop=FALSE]
            rfout$nodepred <-
                rfout$nodepred[1:max.nodes, , drop=FALSE]
            rfout$xbestsplit <-
                rfout$xbestsplit[1:max.nodes, , drop=FALSE]
            rfout$leftDaughter <-
                rfout$leftDaughter[1:max.nodes, , drop=FALSE]
            rfout$rightDaughter <-
                rfout$rightDaughter[1:max.nodes, , drop=FALSE]
        }
        cl <- match.call()
        cl[[1]] <- as.name("rfPoisson")
        ## Make sure those obs. that have not been OOB get NA as prediction.
        ypred <- rfout$ypred
        if (any(rfout$oob.times < 1)) {
            ypred[rfout$oob.times == 0] <- NA
        }
        
        out <- list(call = cl,
                    type = "regression",
                    predicted = structure(ypred, names=x.row.names),
                    dev = rfout$dev,
                    oob.times = rfout$oob.times,
                    importance = if (importance) matrix(rfout$impout, p, 2,
                    dimnames=list(x.col.names,
                                  c("%IncLossFunction","IncNodePurity"))) else
                        matrix(rfout$impout, ncol=1,
                               dimnames=list(x.col.names, "IncNodePurity")),
                    ntree = ntree,
                    mtry = mtry,
                    forest = if (keep.forest)
                    c(rfout[c("ndbigtree", "nodestatus", "leftDaughter",
                              "rightDaughter", "nodepred", "bestvar",
                              "xbestsplit")],
                      list(ncat = ncat), list(nrnodes=max.nodes),
                      list(ntree=ntree), list(xlevels=xlevels)) else NULL,
                    y = y,
                    test = if(testdat) {
                        list(predicted = structure(rfout$ytestpred,
                             names=xts.row.names),
                             dev = if(labelts) rfout$devts else NULL)
                    } else NULL,
                    inbag = if (keep.inbag)
                    matrix(rfout$inbag, nrow(rfout$inbag),
                           dimnames=list(x.row.names, NULL)) else NULL)
    class(out) <- "rfCountData"
    return(out)
    }
