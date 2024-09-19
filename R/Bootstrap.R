#' CI.FDF
#'
#' Computes the Bootstrapped confidence intervals for the PFF model
#'
#' @param fit Ouput fit from \code{FDF} function.
#' @param xarg The \code{CLRdata$x.fine} found in output from \code{PNSDdata_fd} function.
#' @param n.rep  Size of bootstrap samples.
#' @param n.sample Original data sample size (used to limit the range of the bootstrap). 
#'                No need if raw data is large enough. 
#'
#' @return  
#'   A list with: 
#' \itemize{
#'   \item \code{Bootstrap.sample} - an array of dimension length(xarg) x #factors x n.rep 
#'   containing the bootstrapped sample
#'   \item \code{x} - same as \code{xarg} from input
#'   \item \code{lower} - bootstrapped lower CI for all factors in fd format
#'   \item \code{upper} - bootstrapped upper CI for all factors in fd format
#' }
#'    
#' @export
#'
#' @importFrom fda fd
#' @importFrom fda eval.fd
#'
#' @examples
#' library(fda)
#' data(fd.data2011)
#' data(CLRData2011)
#' m<- 2
#' factor_output1 <- FDF(data=fd.data2011,h=24*5,k=m,kmax = 10, 
#'                       kern_type = 'BT')
#' example.bootstrap <- CI.FDF(factor_output1,xarg=CLRData2011$x.fine,
#'                                    n.rep = 10)
#'
#' plot(factor_output1$hat.F, lty = 1)
#' plot(example.bootstrap$lower, lty = 2, add = TRUE)
#' plot(example.bootstrap$upper, lty = 2, add = TRUE)
#'
CI.FDF <- function(fit,
                         xarg,
                         n.sample = NULL,
                         n.rep = 100) {

  # assumes that the input is a factor structured object
  stopifnot("kern_type" %in% names(fit))
  
  ### TODO: MOAR SANITY CHECKS
  
  # lag parameter to be used to estimate long-run cov 
  h <- fit$h 
  if (is.null(n.sample)) {
    n.sample <- max(c(h + 10, 500))
  }
  n.factors <- dim(fit$hat.F$coefs)[2]
  n.basis <- dim(fit$hat.F$coefs)[1]
  residuals.fd <- fit$error.fd
  kern_type <- fit$kern_type

  # allow different values for seasonal data analysis
  sz.fit <- dim(fit$hat.beta)[1]
  sz.fit.res <- dim(residuals.fd$coefs)[2]

  # discrete version of the factors
  Bootstrap.sample <- array(data = NA, dim = c(length(xarg), n.factors, n.rep)) 
  # coefficient array (n.basis x n.factors x n.rep)
  Bootstrap.coeff <- array(data = NA, dim = c(n.basis, n.factors, n.rep))
  
  # time where bootstrap sample starts
  idx.start.bootstrap <- sample(1:(sz.fit - n.sample + 1), n.rep, replace = TRUE) 
  
  # could be refactored to run as a lapply() or something more efficient
  for (j in 1:n.rep) {
    # bootstrap on the index to obtain the residual bootstrap replicate 
    idx.res.bootstrap <- sample(x = 1:sz.fit.res, size = n.sample)
    
    # obtain bootstrap sample 
    coeff.fit.boot <- fit$Xhat.fd$coefs[, idx.start.bootstrap[j]:(idx.start.bootstrap[j] + n.sample - 1)]
    coeff.res.boot <- residuals.fd$coefs[, idx.res.bootstrap]
    
    # coefficients of basis functions of a sample Bootstrap
    coeff.fit.boot <- coeff.fit.boot + coeff.res.boot

    # convert into functional data object
    fit.boot.fd <- fda::fd(coeff.fit.boot, fit$Xhat.fd$basis)

    # apply the FDF function to the replicate
    factor.boot <- FDF(data = fit.boot.fd,
                       h = h,
                       k = n.factors,
                       kern_type = kern_type
                      )

    # fixing sing
    factor.boot$hat.F$coefs <- factor.boot$hat.F$coefs %*% 
      diag(sign(diag(inprod(fit$hat.F, factor.boot$hat.F))))
    Bootstrap.sample[, , j] <- fda::eval.fd(xarg, factor.boot$hat.F)
    Bootstrap.coeff[, , j] <- factor.boot$hat.F$coefs
  }

  sd.coeff <- apply(Bootstrap.coeff, c(1, 2), sd)
  sd.fd <- fda::fd(sd.coeff, fit$Xhat.fd$basis)
  CI.lower <- fda::fd(fit$hat.F$coefs - 1.96 * sd.coeff, fit$Xhat.fd$basis)
  CI.upper <- fda::fd(fit$hat.F$coefs + 1.96 * sd.coeff, fit$Xhat.fd$basis)

  return(
    list(
      Bootstrap.sample = Bootstrap.sample,
      x = xarg,
      lower = CI.lower,
      upper = CI.upper
    )
  )
}

