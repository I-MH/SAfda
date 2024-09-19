#' lrc Long-Run Covariance
#'
#' Estimates the long-run covariance operator, as discussed originally in
#' Aue, A., Rice, G. and Sönmez, O. "Detecting and dating structural breaks in functional data without dimension reduction." Journal of the Royal Statistical Society Series B: Statistical Methodology 80.3 (2018): 509-529.
#'
#' @param datafd Input data, in FD format.
#' @param kern_type Kernel type, switched internally to \code{c("BT", "PR", "SP", "FT", "flat")} which are Bartlett, Parzen, Simple, flat_top or flat kernels.
#' @param h The bandwidth parameter, strictly non-zero.
#'
#' @importFrom fda fd
#' @importFrom fda pca.fd
#' @importFrom fda is.fd
#'
#' @return list
#' @export
#'
#' @examples 
#' rnorm(1)
lrc <- function(datafd, kern_type, h) {

  stopifnot(fda::is.fd(datafd), is.character(kern_type),
            kern_type %in% c("BT", "PR", "SP", "FT", "flat"),
            is.numeric(h), h > 0)
  
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
  CCh0 <- datafd$coefs[, 1:N] %*% t(datafd$coefs[, 1:N])
  CCh <- matrix(0, nbasis, nbasis) + CCh0
  
  for (j in 1:h) {
    AuxCCh <- datafd$coefs[, (1 + j):N] %*% t(datafd$coefs[, 1:(N - j)]) +
      t(datafd$coefs[, (1 + j):N] %*%  t(datafd$coefs[, 1:(N - j)]))
    CCh <- CCh + Kernel(j, h, kerneltype) * AuxCCh
  }
  # for the case where it is not symmetric
  CCh <- (CCh + t(CCh)) / 2 

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

#' lrc_paper Long-Run Covariance, Variant 2
#' 
#' Estimates the long-run covariance, as discussed originally in 
#' [Israel's paper]
#'
#' @param datafd Input data, in FD format.
#' @param kern_type Kernel type, switched internally to \code{c("BT", "PR", "SP", "FT", "flat")} which are Bartlett, Parzen, Simple, flat_top or flat kernels.
#' @param h The bandwidth parameter, strictly non-zero.
#' @param p fill in
#'
#' @return list
#' @export
#'
#' @importFrom fda fd
#' @importFrom fda is.fd
#'
#' @examples rnorm(1)
lrc_paper <- function(datafd, kern_type,h,p) {

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

  Datopca <- fda::pca.fd(datafd, nharm = p)
  lam0 <- Datopca$values
  bb0 <- Datopca$harmonics$coefs   # nbasis x num of factors
  BB <- (1 / lam0[1]) * bb0[, 1] %*% t(bb0[, 1])
  for (l in 2:p) {
    BB <- BB + (1 / lam0[l]) * bb0[, l] %*% t(bb0[, l])
  }
  
  CCh <- matrix(0, nbasis, nbasis)
  for (j in 1:h) {
    AuxCCh <- datafd$coefs[, (1 + j):N] %*%  t(datafd$coefs[, 1:(N - j)]) +
      t(datafd$coefs[, (1 + j):N] %*%  t(datafd$coefs[, 1:(N - j)]))
    CCh <- CCh + Kernel(j, h, kerneltype) * AuxCCh
  }
  CCh <- (CCh + t(CCh)) / 2 # in case it is not symmetric

  xbasis <- datafd$basis
  Jinprod <- inprod(xbasis, xbasis)
  Jaux <- eigen(Jinprod)
  Jhalf <- Jaux$vectors %*% diag(sqrt(Jaux$values)) %*% t(Jaux$vectors)
  Jihalf <- solve(Jhalf)
  Ob <- (1 / N) * Jhalf %*% CCh %*% Jinprod %*% BB %*% Jhalf
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

#' lrc_paper_NoInv Long-Run Covariance, Variant 2, No Inverse
#'
#' Estimates the long-run covariance, as discussed originally in 
#' [Israel's paper]. This is the variant without the inverse.
#'
#' @param datafd Input data, in FD format.
#' @param kern_type Kernel type, switched internally to \code{c("BT", "PR", "SP", "FT", "flat")} which are Bartlett, Parzen, Simple, flat_top or flat kernels.
#' @param h The bandwidth parameter, strictly non-zero.
#' @param replicates Number of replicates (default NULL) for the ...
#'
#' @return list
#' @export
#'
#' @importFrom fda fd
#' @importFrom fda is.fd
#'
#' @examples rnorm(1)
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
