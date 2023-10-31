#' CI_Bootstrap
#'
#' Bootstrapped Confidence Intervals
#'
#' @param data Functional input data.
#' @param xarg The \code{x.fine} CLRdata 
#' @param n.sample Number of samples in the input; used to limit the range of the bootstrap
#' @param n.rep  Number of replicates for the bootstrap
#'
#' @return
#'    matrix: length(xarg) x num_factors containing the sd evaluated at xarg
#'    x: a vector containing values where sd is evaluated
#'    
#' @export
#'
#' @importFrom fda fd
#' @importFrom fda eval.fd
#'
#' @examples
#'  example.bootstrap <-
#'    CI.Bootstrap(data = factor.hat,
#'                 xarg,
#'                 n.sample = 2000,
#'                 n.rep = 100)
#'
#'  plot(factor.hat$hat.F, lty = 1)
#'  plot(example.bootstrap$lower, lty = 2, add = TRUE)
#'  plot(example.bootstrap$upper, lty = 2, add = TRUE)
#'
CI_Bootstrap <- function(data,
                         xarg,
                         n.sample = NULL,
                         n.rep = 100) {

  # assumes that the input is a factor structured object
  stopifnot("kern_type" %in% names(data))
  
  ### TODO: MOAR SANITY CHECKS
  
  # lag parameter to be used to estimate long-run cov 
  h <- data$h 
  if (is.null(n.sample)) {
    n.sample <- max(c(h + 10, 500))
  }
  n.factors <- dim(data$hat.F$coefs)[2]
  n.basis <- dim(data$hat.F$coefs)[1]
  residuals.fd <- data$error.fd
  kern_type <- data$kern_type

  # allow different values for seasonal data analysis
  sz.data <- dim(data$hat.beta)[1]
  sz.data.res <- dim(residuals.fd$coefs)[2]

  # discrete version of the factors
  Bootstrap.sample <- array(data = NA, dim = c(length(xarg), n.factors, n.rep)) 
  # coefficient array (n.basis x n.factors x n.rep)
  Bootstrap.coeff <- array(data = NA, dim = c(n.basis, n.factors, n.rep))
  
  # time where bootstrap sample starts
  idx.start.bootstrap <- sample(1:(sz.data - n.sample + 1), n.rep, replace = TRUE) 
  
  # could be refactored to run as a lapply() or something more efficient
  for (j in 1:n.rep) {
    # bootstrap on the index to obtain the residual bootstrap replicate 
    idx.res.bootstrap <- sample(x = 1:sz.data.res, size = n.sample)
    
    # obtain bootstrap sample 
    coeff.fit.boot <- data$Xhat.fd$coefs[, idx.start.bootstrap[j]:(idx.start.bootstrap[j] + n.sample - 1)]
    coeff.res.boot <- residuals.fd$coefs[, idx.res.bootstrap]
    
    # coefficients of basis functions of a sample Bootstrap
    coeff.data.boot <- coeff.fit.boot + coeff.res.boot

    # convert into functional data object
    data.boot.fd <- fda::fd(coeff.data.boot, data$Xhat.fd$basis)

    # apply the FDF function to the replicate
    factor.boot <- FDF(data = data.boot.fd,
                       h = h,
                       k = n.factors,
                       kern_type = kern_type
                      )

    # fixing sing
    factor.boot$hat.F$coefs <- factor.boot$hat.F$coefs %*% 
      diag(sign(diag(inprod(data$hat.F, factor.boot$hat.F))))
    Bootstrap.sample[, , j] <- fda::eval.fd(xarg, factor.boot$hat.F)
    Bootstrap.coeff[, , j] <- factor.boot$hat.F$coefs
  }

  sd.coeff <- apply(Bootstrap.coeff, c(1, 2), sd)
  sd.fd <- fda::fd(sd.coeff, data$Xhat.fd$basis)
  CI.lower <- fda::fd(data$hat.F$coefs - 1.96 * sd.coeff, data$Xhat.fd$basis)
  CI.upper <- fda::fd(data$hat.F$coefs + 1.96 * sd.coeff, data$Xhat.fd$basis)

  return(
    list(
      Bootstrap.sample = Bootstrap.sample,
      x = xarg,
      lower = CI.lower,
      upper = CI.upper
    )
  )
}

