#' Estimation of the Lambda operator - Long-run covariance
#' 
#' Estimates the long-run covariance operator lambda described in the main paper. 
#' 
#' This uses the method discussed originally in
#' Rice, G., and Shang, H.L. (2017), "A Plug-in Bandwidth Selection Procedure for Long-run 
#' Covariance Estimation with Stationary Functional Time Series", 
#' Journal of Time Series Analysis, 38, 591–609.
#'
#'
#' @param datafd Input data, in FD format.
#' @param kern_type Kernel type, switched internally to \code{c("BT", "PR", "SP", "FT", "flat")} which are Bartlett, Parzen, Simple, flat_top or flat kernels.
#' @param h The bandwidth parameter, strictly non-zero.
#' @param replicates Number of replicates (default NULL) for the ...
#'
#' @return 
#' A list with:
#' \itemize{
#'   \item \code{eigenf} - eigenfunctions of the operator
#'   \item \code{eigenval} - eigenvalues of the operator
#'   \item \code{lrc.k} - kernel of the lambda operator
#' }
#' @export
#'
#' @importFrom fda fd
#' @importFrom fda is.fd
lrc_paper_NoInv <- function(datafd, 
                            kern_type,
                            h, 
                            replicates = NULL) {

  stopifnot(fda::is.fd(datafd), is.character(kern_type),
            kern_type %in% c("BT", "PR", "SP", "FT", "flat"),
            is.numeric(h), h > 0)
 
  # returns the eigenfunctions of the long run cov
  kerneltype <- switch(
    kern_type,
    BT = "Bartlett",
    PR = "Parzen",
    FT = "flat_top",
    SP = "Simple",
    flat = "flat"
  )

  N <- dim(datafd$coefs)[2]
  nbasis <- datafd$basis$nbasis
  if (is.null(replicates)) {
    CCh <- matrix(0, nbasis, nbasis)
    for (j in 1:h) {
      AuxCCh <- datafd$coefs[, (1 + j):N] %*% t(datafd$coefs[, 1:(N - j)]) +
        t(datafd$coefs[, (1 + j):N] %*%  t(datafd$coefs[, 1:(N - j)]))
      CCh <- CCh + Kernel(j, h, kerneltype) * AuxCCh
    }
    # for the case that it is not symmetric
    CCh <- (CCh + t(CCh)) / 2 
  } else {
    nrep <- length(unique(replicates))
    CChrep <- array(0, c(nbasis, nbasis, nrep))
    AllCoeff <- datafd$coefs
    for (k in 1:nrep) {
      Coeff <- AllCoeff[, replicates == k]
      N2 <- dim(Coeff)[2]
      for (j in 1:h) {
        AuxCCh <- Coeff[, (1 + j):N2] %*% t(Coeff[, 1:(N2 - j)]) +
          t(Coeff[, (1 + j):N2] %*%  t(Coeff[, 1:(N2 - j)]))
        CChrep[, , k] <- CChrep[, , k] + Kernel(j, h, kerneltype) * AuxCCh
      }
      # for the case that it is not symmetric
      CChrep[, , k] <- (CChrep[, , k] + t(CChrep[, , k])) / 2
    }
    CCh <-  apply(CChrep, 1:2, mean)
    # and one final non-symmetric catch
    CCh <- (CCh + t(CCh)) / 2
  }

  xbasis <- datafd$basis
  Jinprod <- inprod(xbasis, xbasis)
  Jaux <- eigen(Jinprod)
  Jhalf <- Jaux$vectors %*% diag(sqrt(Jaux$values)) %*% t(Jaux$vectors)
  Jihalf <- solve(Jhalf)
  Ob <- (1 / N) * Jhalf %*% CCh %*% Jhalf
  result <- eigen(Ob)
  
  # coefficients of the eigenfunctions
  bb <- Jihalf %*% Re(result$vectors) 
  # eigenvalues of the operator
  eigenval <- Re(result$values) 
  eigenf <- fda::fd(Re(bb), basisobj = xbasis)
  
  return(list(
    lrc.k = CCh,
    Ob = Ob,
    bb = Re(bb),
    eigenf = eigenf,
    eigenval = eigenval
  ))
}
