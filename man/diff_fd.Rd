% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{diff_fd}
\alias{diff_fd}
\title{diff_fd}
\usage{
diff_fd(datafd)
}
\arguments{
\item{datafd}{Input data in FD format.}
}
\value{
Differentiated (1-step differencing operator) input, in `fd` format.
}
\description{
1-step difference of a data set in class `fd` format, returned as `fd`.
}
\examples{
rnorm(1)
}
