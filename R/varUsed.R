#' varUsed
#' @name varUsed
#' @description Find out which predictor variables are actually used in the random forest. 
#' @param x, An object of class randomForest.
#' @param by.tree, Should the list of variables used be broken down by trees in the forest?
#' @param count, Should the frequencies that variables appear in trees be returned?
#' @return If \code{count=TRUE} and \code{by.tree=FALSE}, a integer vector containing frequencies 
#' that variables are used in the forest. If \code{by.tree=TRUE}, a matrix is returned, breaking down 
#' the counts by tree (each column corresponding to one tree and each row to a variable).\cr \cr
#' If \code{count=FALSE} and \code{by.tree=TRUE}, a list of integer indices is returned giving 
#' the variables used in the trees, else if \code{by.tree=FALSE}, a vector of integer indices 
#' giving the variables used in the entire forest. 
#' @author Andy Liaw
#' @export
varUsed <- function(x, by.tree=FALSE, count=TRUE) {
    if (!inherits(x, "rfCountData"))
        stop(deparse(substitute(x)), "is not a rfCountData object")
    if (is.null(x$forest))
        stop(deparse(substitute(x)), "does not contain forest")
    
    p <- length(x$forest$ncat)  # Total number of variables.
    if (count) {
        if (by.tree) {
            v <- apply(x$forest$bestvar, 2, function(x) {
                xx <- numeric(p)
                y <- table(x[x>0])
                xx[as.numeric(names(y))] <- y
                xx
            })
        } else {
            v <- numeric(p)
            vv <- table(x$forest$bestvar[x$forest$bestvar > 0])
            v[as.numeric(names(vv))] <- vv
        }
    } else {
        v <- apply(x$forest$bestvar, 2, function(x) sort(unique(x[x>0])))
        if(!by.tree) v <- sort(unique(unlist(v)))
    }
    v
}
