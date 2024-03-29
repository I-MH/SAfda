% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FDF.R
\name{FDF}
\alias{FDF}
\title{FDF}
\usage{
FDF(
  argvals = NULL,
  data,
  stationary = TRUE,
  h = 5,
  k = NULL,
  kmax = 6,
  p = 5,
  kern_type = "BT",
  rangeval = c(0, 1),
  nbasis = 25,
  basis = "Fourier",
  lambda = NULL,
  replicates = NULL,
  plot = FALSE
)
}
\arguments{
\item{argvals}{A vector containing the points where the function values are observed. If NULL, assumed to be equally spaced in rangeval.}

\item{data}{A matrix, `m x N`, where `N` is the number of curves.}

\item{stationary}{Logical, TRUE implies stationarity.}

\item{h}{Number of lags to be used to compute the long-run covariance.}

\item{k}{Number of factors to be used. If NULL, then estimated as in the reference.}

\item{kmax}{Maximum number of factors to be tested.}

\item{p}{Number of eigenfunctions to be used in approximating the inverse covariance operator.}

\item{kern_type}{Type of kernel to be used when estimating the long run covariance operator.}

\item{rangeval}{Numeric, length 2, defining the interval over which the functional data object can be evaluated.}

\item{nbasis}{Number of basis functions to be used for the functional data.}

\item{basis}{Type of basis, from options 'Fourier' and 'Bspline'}

\item{lambda}{Non-negative real number, smoothing parameter to be applied to the estimated functional parameter. If NULL, estimated using generalized CV.}

\item{replicates}{???}

\item{plot}{Logical, plot to be displayed or not.}
}
\value{
list()
  hat.beta: the estimated time series
  hat.F: the estimated F
  hat.K_ratio: the estimated number of factors
  Xhat: fitted values of the data
  eigenval: eigenvalues of the long run cov
  eigenvalC0: eigenvalues of the cov at lag 0
  Ob: (not for users)
  Ob_result: (not for users)
}
\description{
FDF
}
\examples{
rnorm(1)
}
