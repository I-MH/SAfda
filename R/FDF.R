#' FDF
#'
#' Fits a functional dynamic factor model for Hilbert space data, can be applied to densities data after 
#' preprocessing with PNSDdata_fd
#' 
#' @param argvals  A vector containing the points where the function values are observed. If NULL, assumed to be equally spaced in rangeval.
#' @param data A matrix, `m x N`, where `N` is the number of curves.
#' @param stationary Logical, TRUE implies stationarity.
#' @param h  Number of lags to be used to compute the long-run covariance.
#' @param k  Number of factors to be used. If NULL, then estimated as in the reference.
#' @param kmax  Maximum number of factors to be tested.
#' @param p  Number of eigenfunctions to be used in approximating the inverse covariance operator.
#' @param kern_type  Type of kernel to be used when estimating the long run covariance operator.
#' @param rangeval  Numeric, length 2, defining the interval over which the functional data object can be evaluated.
#' @param nbasis Number of basis functions to be used for the functional data.
#' @param basis Type of basis, from options 'Fourier' and 'Bspline'
#' @param lambda  Non-negative real number, smoothing parameter to be applied to the estimated functional parameter. If NULL, estimated using generalized CV.
#' @param replicates  ???
#' @param plot  Logical, plot to be displayed or not.
#'
#' @importFrom fda fd
#'
#' @return  list()
#'   hat.beta: the estimated time series
#'   hat.F: the estimated F
#'   hat.K_ratio: the estimated number of factors
#'   Xhat: fitted values of the data
#'   eigenval: eigenvalues of the long run cov
#'   eigenvalC0: eigenvalues of the cov at lag 0
#'   Ob: (not for users)
#'   Ob_result: (not for users)
#' @export
#'
#' @examples 
#' data(fd.data2011)
#' m<- 5
#' factor_output1 <- FDF(data=fd.data2011,h=24*5,k=m,kmax = 10, kern_type = 'BT',plot=TRUE)
FDF <- function(argvals = NULL,
                data,
                stationary = TRUE,
                h = 5,
                k = NULL,
                kmax = 6,
                p = 5,
                kern_type = "BT",
                rangeval = c(0, 1),
                nbasis = 25,
                basis = 'Fourier',
                lambda = NULL,
                replicates = NULL,
                plot = FALSE) {


  # transform data to fd class ------------------------------
  newData <- data.fd(argvals = argvals,
                     data = data,
                     basis = basis,
                     nbasis = nbasis,
                     rangeval = rangeval,
                     lambda = lambda
                    )

  datafd0 <- newData$datafd
  rangeval <- newData$rangeval
  N <- newData$N
  m <- newData$m
  nbasis <- newData$nbasis
  CoeffData <- t(coef(datafd0))    # N x nbasis
  xbasis <- datafd0$basis
  Jinprod <- inprod(xbasis, xbasis)

  if (stationary) {
    datafd <- datafd0
  } else {
    datafd <- diff_fd(datafd = datafd0)
  }

  # estimate long run cov -----------------------------------
  result.lrc <- lrc_paper_NoInv(datafd = datafd,
                                kern_type = kern_type,
                                h = h,
                                replicates = replicates
                               )

  # estimate number of factors
  hatK_aux <- hatK(Re(result.lrc$eigenval[1:kmax]))
  if (is.null(k)) {
    k <- hatK_aux$hat_k
  }
  rank_bb <- qr(result.lrc$Ob)$rank
  newr <- min(k, rank_bb)

  # estimating factors
  bb.f <- result.lrc$bb[, 1:newr] # coef for factor
  e.val <- Re(result.lrc$eigenval[1:newr]) #eigenvalues of the operator
  #loading factors
  lam <- fda::fd(Re(bb.f), basisobj = xbasis)
  #factor processes
  ff <- CoeffData %*% Jinprod %*% bb.f

  # fitting
  if (is.null(argvals)) {
    tt <- seq(rangeval[1], rangeval[2], length.out = m)
  } else {
    tt = argvals
  }
  Xhat <-  eval.fd(tt, lam) %*% t(as.matrix(ff))

  # fd version
  coeff.fit <- lam$coefs %*% t(ff)
  Xhat.fd <- fda::fd(coeff.fit, basisobj = xbasis)
  error.fd <- fda::fd(datafd0$coefs - coeff.fit, basisobj = xbasis)

  # obtain scree plot if desired
  hatk0 <- which.min(Re(result.lrc$eigenval[1:kmax]))
  if(plot) {
    par(mfrow = c(2, 1))
    plot(Re(result.lrc$eigenval[1:kmax]),
         type = 'o',
         main = 'Scree Plot',
         ylab = 'eigenvalues',
         xlab = 'Number of Sources'
       )
    plot(hatK_aux$ratio,
         type = 'o',
         main = 'Ratio Scree Plot',
         ylab = 'ratio',
         xlab = 'Number of Sources'
        )
    par(mfrow = c(1, 1))
  }

  # make lam.fd a pca.fd class for varmx.pca.fd()
  lam.fd <- list(lam, result.lrc$eigenval, ff)
  class(lam.fd) <- "pca.fd"
  names(lam.fd) <- c("harmonics", "values", "scores")
  ff.Varmx <- varmx.pca.fd(lam.fd, nharm = dim(ff)[2])

  return(
    list(
      hat.beta = ff,
      hat.F = lam,
      Varmx = ff.Varmx,
      hat.beta.Varmx = ff.Varmx$scores,
      hat.f.Varmx = ff.Varmx$harmonics,
      hat.K_ratio = hatK_aux$hat_k,
      hat.K.scree = hatk0,
      h = h,
      kern_type = kern_type,
      Xhat = Xhat,
      Xhat.fd = Xhat.fd,
      error.fd = error.fd,
      eigenval = Re(result.lrc$eigenval),
      argvals = tt
    )
  )
}
