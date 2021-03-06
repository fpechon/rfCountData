% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importance.R
\name{importance}
\alias{importance}
\title{Extract variable importance measure}
\usage{
importance(x)
}
\arguments{
\item{x, }{an object of class \code{\link{rfCountData}}.}
}
\value{
A vector of importance measure, one item for each predictor variable.
}
\description{
Extract variable importance measure
}
\details{
The measure is computed from permuting OOB data: For each tree, the prediction error 
on the out-of-bag portion of the data is recorded (mean deviance). 
Then the same is done after permuting each predictor variable. 
The difference between the two are then averaged over all trees. Note that 
the mean deviances are compared (and not the deviances themselves).
}
\seealso{
\link{rfPoisson}, \link{varImpPlot}.
}
