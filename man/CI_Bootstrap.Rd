% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bootstrap.R
\name{CI_Bootstrap}
\alias{CI_Bootstrap}
\title{CI_Bootstrap}
\usage{
CI_Bootstrap(data, xarg, n.sample = NULL, n.rep = 100)
}
\arguments{
\item{data}{Functional input data.}

\item{xarg}{The \code{x.fine} CLRdata}

\item{n.sample}{Number of samples in the input; used to limit the range of the bootstrap}

\item{n.rep}{Number of replicates for the bootstrap}
}
\value{
matrix: length(xarg) x num_factors containing the sd evaluated at xarg
   x: a vector containing values where sd is evaluated
}
\description{
Bootstrapped Confidence Intervals
}
\examples{
 example.bootstrap <-
   CI.Bootstrap(data = factor.hat,
                xarg,
                n.sample = 2000,
                n.rep = 100)

 plot(factor.hat$hat.F, lty = 1)
 plot(example.bootstrap$lower, lty = 2, add = TRUE)
 plot(example.bootstrap$upper, lty = 2, add = TRUE)

}
