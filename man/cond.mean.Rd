% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/condmoments.R
\name{cond.mean}
\alias{cond.mean}
\title{Conditional mean estimation by kernel smoothing}
\usage{
cond.mean(datax, datay, x0, h)
}
\arguments{
\item{datax}{vector or matrix of observed X values (conditioning variables).}

\item{datay}{vector or matrix of observed Y values (conditioned variables).
Must have the same number of observations as \code{datax}}

\item{x0}{point at which the U-statistic should be computed.}

\item{h}{the bandwidth of the kernel smoothing.}
}
\value{
a vector having \code{NROW(datay)} entries with the conditional mean
of the corresponding variables.
}
\description{
Conditional mean estimation by kernel smoothing
}
\examples{
n = 500
datax = rnorm(n)
datay = datax + rnorm(n)
x0 = 1
h = 0.2
cond.mean(datax = datax, datay = datay, x0 = x0, h = h)

}
