% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rfPoisson.R
\name{rfPoisson}
\alias{rfPoisson}
\title{rfPoisson}
\usage{
rfPoisson(x, offset = NULL, y = NULL, xtest = NULL, ytest = NULL,
  offsettest = NULL, ntree = 500, mtry = max(floor(ncol(x)/3), 1),
  replace = TRUE, sampsize = if (replace) nrow(x) else ceiling(0.632 *
  nrow(x)), nodesize = 5000, maxnodes = NULL, importance = TRUE,
  nPerm = 1, do.trace = FALSE, keep.forest = !is.null(y) &&
  is.null(xtest), keep.inbag = FALSE, ...)
}
\arguments{
\item{x}{a data frame or a matrix of predictors.}

\item{offset}{a vector of same size as y, corresponding to the log of observation time (e.g. log of exposure). Default is 0.}

\item{y}{a vector of Poisson responses}

\item{xtest}{a data frame or matrix (like \code{x}) containing predictors for the test set.}

\item{ytest}{response for the test set}

\item{offsettest}{Offset for the test set (like \code{offset})}

\item{ntree, }{Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.}

\item{mtry}{Number of variables randomly sampled as candidates at each split.  Default is p/3, where p is the number of variables in \code{x}.}

\item{replace}{Should sampling of cases be done with or without replacement?}

\item{sampsize}{Size(s) of sample to draw.}

\item{nodesize}{Minimum size of terminal nodes.  Setting this number 
larger causes smaller trees to be grown (and thus take less time). Default is 5000.}

\item{maxnodes}{Maximum number of terminal nodes trees in the forest can have.
If not given, trees are grown to the maximum possible 
(subject to limits by \code{nodesize}).  If set larger than maximum possible, a warning is issued.}

\item{importance}{Should importance of predictors be assessed? Default is \code{TRUE}.}

\item{nPerm}{Not yet implemented. Number of times the OOB data are permuted per tree for 
assessing variable importance.  Number larger than 1 gives slightly more stable estimate, 
but not very effective.}

\item{do.trace}{If set to \code{TRUE}, give a more verbose output as \code{rfPoisson} is run.
If set to some integer, then running output is printed for every \code{do.trace} trees.}

\item{keep.forest}{If set to \code{FALSE}, the forest will not be retained in the output object.
If \code{xtest} is given, defaults to \code{FALSE}.}

\item{keep.inbag}{Should an \code{n} by \code{ntree} matrix be returned that keeps track of which 
samples are `in-bag' in which trees (but not how many times, if sampling with replacement}

\item{...}{other parameters passed to lower functions.}
}
\value{
TBC
}
\description{
\code{rfPoisson} implements Breiman's random forest algorithm (based on 
Breiman and Cutler's original Fortran code) for regression and has been modified 
to be used with Poisson data that have different observation periods.\cr
More specifically, the best split is the one that will maximise the decrease of the poisson deviance. An offset
has also been introduced to accomodate for different times of exposure. The offset the log of the exposure.
}
\note{
TBC
}
\examples{
if (!require(CASdatasets)) install.packages("CASdatasets", repos = "http://cas.uqam.ca/pub/R/", type="source")
require(CASdatasets)
data("freMTPLfreq")
library(rfCountData)
m0 = rfPoisson(y = freMTPLfreq[1:10000,]$ClaimNb,
                  offset = log(freMTPLfreq[1:10000,]$Exposure),
                  x = freMTPLfreq[1:10000,c("Region", "Power", "DriverAge")],
                  ntree = 20)
predict(m0, newdata = freMTPLfreq[10001:10050,c("Region", "Power", "DriverAge")], 
offset = log(freMTPLfreq[10001:10050,"Exposure"]))
}
\references{
R package \emph{randomForest},  \url{https://cran.r-project.org/package=randomForest}\cr
Breiman, L. (2001), \emph{Random Forests}, Machine Learning 45(1), 5-32.\cr
Breiman, L (2002), ``Manual On Setting Up, Using, And Understanding Random Forests V3.1'', \url{https://www.stat.berkeley.edu/~breiman/Using_random_forests_V3.1.pdf}.
}
\author{
Florian Pechon, \email{florian.pechon@uclouvain.be}, based on the package randomForest by
 Andy Liaw \email{andy_liaw@merck.com} and Matthew Wiener \email{matthew_wiener@merck.com} based 
 on original Fortran code by Leo Breiman and Adele Cutler.
}
