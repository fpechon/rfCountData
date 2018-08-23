#' treesize
#' @param x, an object of class \code{rfCountData}, which contains a \code{forest} component.
#' @param terminal, count terminal nodes only (\code{TRUE}) or all nodes (\code{FALSE}
#' @return A vector containing number of nodes for the trees in the \code{rfCountData} object.
#' @note The \code{rfCountData} object must contain the \code{forest} component; i.e., created with \code{rfPoisson(...,keep.forest=TRUE)}.
#' @author Andy Liaw \email{andy_liaw@merck.com}
#' @export
treesize <- function(x, terminal=TRUE) {
  if(!inherits(x, "rfCountData"))
    stop("This function only works for objects of class `rfCountData'")
  if(is.null(x$forest)) stop("The object must contain the forest component")
  if(terminal) return((x$forest$ndbigtree+1)/2) else return(x$forest$ndbigtree)
}
