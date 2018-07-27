#' getTree
#' @description This function extract the structure of a tree from a \code{rfCountData} object.
#' @usage getTree(rfobj, k=1, labelVar=FALSE)
#' @param rfobj, a \code{rfCountData} object.
#' @param k, which tree to extract ?
#' @param labelVar, Should better labels be used for splitting variables and predicted class ?
#' @details For numerical predictors, data with values of the variable less than or equal to the 
#' splitting point go to the left daughter node. \cr
#' For categorical predictors, the splitting point is represented by an integer, whose binary 
#' expansion gives the identities of the categories that goes to left or right. \cr
#' For example, if a predictor has four categories, and the split point is 13. 
#' The binary expansion of 13 is (1, 0, 1, 1) (because $13 = 1*2^0 + 0*2^1$ + $1*2^2 + 1*2^3$),
#' so cases with categories 1, 3, or 4 in this predictor get sent to the left, and the rest to the right.
#' @return A matrix (or data frame, if \code{labelVar=TRUE}) with six columns 
#' and number of rows equal to total number of nodes in the tree.  The six columns are: \cr
#' \code{left daughter}, the row where the left daughter node is; 0 if the node is terminal.\cr
#' \code{right daughter}, the row where the right daughter node is; 0 if the node is terminal.\cr
#' \code{split var}, which variable was used to split the node; 0 if the  node is terminal.\cr
#' \code{split point}, where the best split is; see Details for  categorical predictor.\cr
#' \code{status}, is the node terminal (-1) or not (1).\cr
#' \code{prediction},the prediction for the node; 0 if the node is not terminal.
#' @author Andy Liaw \email{andy_liaw@merck.com}.
#' @seealso \link{rfPoisson}
#' @export
getTree <- function(rfobj, k=1, labelVar=FALSE) {
  if (is.null(rfobj$forest)) {
    stop("No forest component in ", deparse(substitute(rfobj)))
  }
  if (k > rfobj$ntree) {
    stop("There are fewer than ", k, "trees in the forest")
  }
  if (rfobj$type == "regression") {
      tree <- cbind(rfobj$forest$leftDaughter[,k],
                    rfobj$forest$rightDaughter[,k],
                    rfobj$forest$bestvar[,k],
                    rfobj$forest$xbestsplit[,k],
                    rfobj$forest$nodestatus[,k],
                    rfobj$forest$nodepred[,k])[1:rfobj$forest$ndbigtree[k],]
  } else {
      tree <- cbind(rfobj$forest$treemap[,,k],
                    rfobj$forest$bestvar[,k],
                    rfobj$forest$xbestsplit[,k],
                    rfobj$forest$nodestatus[,k],
                    rfobj$forest$nodepred[,k])[1:rfobj$forest$ndbigtree[k],]
  }

  dimnames(tree) <- list(1:nrow(tree), c("left daughter", "right daughter",
                                         "split var", "split point",
                                         "status", "prediction"))

  if (labelVar) {
      tree <- as.data.frame(tree)
      v <- tree[[3]]
      v[v == 0] <- NA
      tree[[3]] <- factor(rownames(rfobj$importance)[v])
      if (rfobj$type == "classification") {
          v <- tree[[6]]
          v[! v %in% 1:nlevels(rfobj$y)] <- NA
          tree[[6]] <- levels(rfobj$y)[v]
      }
  }
  tree
}

