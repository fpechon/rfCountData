% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partialPlot.R
\name{partialPlot}
\alias{partialPlot}
\title{Partial dependence plot}
\usage{
partialPlot(x, pred.data, x.var, offset,
w, plot = TRUE, n.pt = min(length(unique(pred.data[, xname])), 51),
rug = TRUE, xlab=deparse(substitute(x.var)), ylab="",
main=paste("Partial Dependence on", deparse(substitute(x.var))),...)
}
\arguments{
\item{x, }{an object of class \code{rfCountData}, which contains a \code{forest} component.}

\item{pred.data, }{a data frame used for contructing the plot, usually the training data used to contruct the random forest.}

\item{x.var, }{name of the variable for which partial dependence is to be examined.}

\item{offset, }{a vector, the corresponding log-exposures of pred.data.}

\item{w, }{weights to be used in averaging; if not supplied, mean is not weighted}

\item{plot, }{whether the plot should be shown on the graphic device.}

\item{n.pt, }{if \code{x.var} is continuous, the number of points on the grid for evaluating partial dependence.}

\item{rug, }{whether to draw hash marks at the bottom of the plot indicating the deciles of \code{x.var}.}

\item{xlab, }{label for the x-axis.}

\item{ylab, }{label for the y-axis.}

\item{main, }{main title for the plot.}

\item{..., }{other graphical parameters to be passed on to \code{plot} or \code{lines}.}
}
\value{
A list with two components: \code{x} and \code{y}, which are the values used in the plot. 
The \code{rfCountData} object must contain the forest component; i.e., 
created with \code{rfPoisson(..., keep.forest=TRUE)}. This function runs quite slow for large data sets.
}
\description{
Partial dependence plot gives a graphical
 depiction of the marginal effect of a variable on the response variable.
}
\details{
The function being plotted is defined as: 
\deqn{f(x) = 1/n* \sum_{i=1}^{n} f(x, x_{iC})},
where x is the variable for which partial dependence is sought, and x_{iC} is the other variables in the data.
}
\references{
Friedman, J. (2001). Greedy function approximation: the gradient boosting machine, Ann. of Stat.
}
\seealso{
\link{rfPoisson}
}
\author{
Andy Liaw \email{andy_liaw@merck.com}
}
