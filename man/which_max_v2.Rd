% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{which_max_v2}
\alias{which_max_v2}
\title{which_max_v2 An updated `which.max`}
\usage{
which_max_v2(z)
}
\arguments{
\item{z}{Input numeric vector.}
}
\value{
a `which.max` one-in inset from the boundary.
}
\description{
Computes `which.max`, but remaps boundary values to
one-in inset values.
}
\examples{
x <- sample(1:10, 10, TRUE)
which_max_v2(x)

}
