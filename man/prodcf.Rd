% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{prodcf}
\alias{prodcf}
\title{prodcf}
\usage{
prodcf(c, f, s)
}
\arguments{
\item{c}{Power}

\item{f}{Base}

\item{s}{Observation times; first difference is the base for the sum.}
}
\value{
A list of two objects: the `sum` and the `norm` (constant)
}
\description{
Utility function to compute the midpoint-rule Riemann sum of `f^c`
with possibly irregular bases (from the first-difference of `s`) and 
then return the product of `f^c`, normalized by that sum.
#'
}
\examples{
rnorm(1)
}
