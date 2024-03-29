% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hatK.R
\name{hatK}
\alias{hatK}
\title{hatK}
\usage{
hatK(evalues)
}
\arguments{
\item{evalues}{Numeric vector of eigenvalues}
}
\value{
list
hat_k  Index of the best k
ratio_k  The value of the best k
}
\description{
Computes the auxiliary K, as the ratio of each eigenvalue to its predecessor;
returns the index and value of the minimum of such quotients.
}
\examples{
rnorm(1)
}
