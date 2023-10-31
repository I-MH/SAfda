#
#  this file is largely done and refactored; just need to check the docs, and write unit tests,
#  and put in sanity checks at the start
#
#' z_spline_basis: Computes a basis set of Z-splines
#' @description Method described in Hron, Menafoglio, Templ, Hruzova and Filzmoser (2016),
#' `Simplicial principal component analysis for density functions in Bayes spaces`,
#' Computational Statistics & Data Analysis, Volume 94, Pages 330-350, ISSN 0167-9473.
#' @param knots Input knot list for spline basis creation.
#' @param k Order of the B-spline basis, one higher than their degree.
#' @importFrom fda create.bspline.basis
#' @importFrom fda eval.basis
#' @return list(C0 = C0, M0 = M, K = L, D = Dmat)
#'   Four matrix objects.
#' @export
#'
#' @examples
#'   knots <- seq(from = 1, to = 20, length.out = 100)
#'   k <- 4
#'   res <- z_spline_basis(knots = knots, k = k)
z_spline_basis <- function(knots, k) {
  stopifnot(length(knots) > 2)

  r <- length(knots)
  lambda_index <- 0:(r-1)
  g <- lambda_index[length(lambda_index) - 1]
  N <- g + (k-1) + 1

  lambda <- c(rep(min(knots), k-1),
              knots,
              rep(max(knots), k-1))
  div <- seq(from = min(lambda), to = max(lambda), length.out = 1000)

  # standard B-spline basis; collocation matrix := C
  spline_basis <- fda::create.bspline.basis(rangeval = range(knots),
                                            nbasis = N,
                                            norder = k,
                                            breaks = knots)
  Cmat <- fda::eval.basis(evalarg = div,
                          basisobj = spline_basis)

  # Matrix D
  differ <- lambda[(1+k):(r + 2 * (k-1))] - lambda[1:(r + k -2)]
  Dmat <- k * diag(1 / differ)

  # Matrix L
  L <- array(data = 0, dim = c(N, N-1))
  diag(L) <- 1 # diagonal
  diag(L[-1, ]) <- -1  # first subdiagonal

  # Spline0 basis: collocation matrix C0
  C0 <- Cmat %*% Dmat %*% L

  # Matrix M
  division <- seq(from = min(lambda), to = max(lambda), length.out = 10000)
  step <- diff(div[1:2])
  CC <- eval.basis(division, spline_basis)

  CC0 <- CC %*% Dmat %*% L

  # refactored: filling the upper triangle is still clumsy
  M <- diag(rep(0, N - 1))
  slp <- c(0.5, rep(1, nrow(CC0) - 2), 0.5)
  for(i in 1:(N-1)) {
    for(j in i:(N-1)) { # only do the upper triangle
      M[i,j] <- step * (CC0[, i] * CC0[, j]) %*% slp
    }
  }
  M[lower.tri(M, diag = FALSE)] <- t(M)[lower.tri(M, diag = FALSE)]

  return(list(C0 = C0,
              M0 = M,
              K = L,
              D = Dmat))
}
