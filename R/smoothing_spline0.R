#
#  this file is largely done and refactored; just need to check the docs, and write unit tests,
#  and put in sanity checks at the start
#
#' @title Modified Smoothing Spline
#'
#' @description Method described in Hron, Menafoglio, Templ, Hruzova and Filzmoser (2016),
#' `Simplicial principal component analysis for density functions in Bayes spaces`,
#' Computational Statistics & Data Analysis, Volume 94, Pages 330-350, ISSN 0167-9473.'
#' Extension of the original grid by adding (k-1) knots.
#'
#' @param knots Input knot list for spline basis creation.
#' @param tp Points of approximation.
#' @param f Values at tp.
#' @param w Weight coefficients.
#' @param k Order of the B-spline basis, one higher than their degree.
#' @param der Derivation?
#' @param alpha Smoothing parameter.
#' @param ch if FALSE, then functional with (1-alpha) and alpha; if TRUE, then functional with alpha alone.
#' @importFrom splines splineDesign
#'
#' @return A list of three objects: `list(J, z, spline0))`
smoothing_spline0 = function(knots,
                            tp,
                            f,
                            w,
                            k,
                            der,
                            alpha,
                            ch = FALSE) {

  r <- length(knots)

  # add copies of first and last knots out k-1 in both directions
  y <- c(rep(knots[1], k-1),
         knots,
         rep(knots[length(knots)], k-1))

  # Collocation matrix K
  K <- splines::splineDesign(y, tp, k, outer.ok = TRUE)

  # Diag matrix with weights
  W <- diag(w)

  # Collocation matrix C
  div <- seq(from = min(y), to = max(y), length.out = 1000)
  lambda <- c(0:(r - 1))
  g <- lambda[length(lambda) - 1]

  # Dimension (space of splines)
  N <- g + (k - 1) + 1
  Cmat <- array(data = 0, dim = c(length(div), N))
  for(i in 1:N) {
    Cmat[, i] <- splineDesign(knots = y[i - 1 + 1:(k + 1)],
                              x = div, ord = k, outer.ok = TRUE)
  }

  # Verification of full column rank of collocation matrix K
  if (length(tp) <= N) stop('length(t) must be higher then Dimension(space of splines)')
  if (qr(K)$rank != N) stop('Collocation matrix does not have full column rank.')

  # Matrix S
  if(der == 0) {
    S <- diag(x = 1, nrow = N, ncol = N)
  } else if(der > 0) {
    i <- der
    while (i > 0) {
      rozdil <- y[(1+k):(N+k-i)] - y[(1+i):(N)]
      Dmat <- (k - i) * diag(1 / rozdil)
      Lmat <- array(data = 0, dim = c(N - i, N - i + 1))
      diag(Lmat) <- -1 # diagonal
      diag(Lmat[, -1]) <- 1  # first superdiagonal
      if (i == der) { # initialize
        S = Dmat %*% Lmat
      } else { # composition
        S = S %*% Dmat %*% Lmat
      }
      i=i-1
    } # iterate downward
  }

  # Matrix M - augment knot sequence - order of spline = k - der
  kk = k - der
  Y <- c(rep(knots[1], kk - 1),
         knots,
         rep(knots[length(knots)], kk - 1))

  division <- seq(from = min(Y), to = max(Y), length.out = 10000)
  Lambda <- 0:(r-1)
  G <- Lambda[r - 1]

  # Matrix M - spline space dimension
  NN <- G + (kk - 1) + 1

  # Matrix M - collocation matrix KK
  CC <- splineDesign(knots = Y, x = division, ord = kk, outer.ok = TRUE)

  step <- diff(division[1:2])

  # Matrix M
  M <- diag(rep(0, NN))
  slp <- c(0.5, rep(1, nrow(CC) - 2), 0.5)
  for(i in 1:NN) {
    for(j in i:NN) { # only do the upper triangle
      M[i,j] <- step * (CC[, i] * CC[, j]) %*% slp
    }
  }
  # fill in the lower.tri by transposing then copying
  M[lower.tri(M, diag = FALSE)] <- t(M)[lower.tri(M, diag = FALSE)]

  # Matrix D
  differ <- y[(1 + k):(r + 2 * (k - 1))] - y[(1:(r + k - 2))]
  Dmat <- k * diag(1 / differ)

  # Matrix K
  KK <- array(data = 0, dim = c(N, N - 1))
  diag(KK) <- 1 # diagonal
  diag(KK[-1, ]) <- -1  # first subdiagonal

  # Matrix U
  U = S %*% Dmat %*% KK

  # Matrix G, vector g
  GG <- alpha * t(KK) %*% t(Dmat) %*% t(K) %*% W %*% K %*% Dmat %*% KK +
  if(!ch) {
    (1 - alpha) * t(U) %*% M %*% U
  } else {
                 t(U) %*% M %*% U
  }

  gg <- alpha * t(KK) %*% t(Dmat) %*% t(K) %*% W %*% f

  # vector of B-spline coefficients := z
  z <- solve(GG) %*% gg

  # B-spline basis in L20
  Bbasis <- Cmat %*% Dmat %*% KK

  # Resulting spline
  spline0 <- (Cmat %*% Dmat %*% KK) %*% z
  integral <- step * t(spline0) %*% c(0.5, rep(1, length(spline0) - 2), 0.5)

  J <- alpha * t(f - K %*% Dmat %*% KK %*% z) %*% W %*% (f - K %*% Dmat %*% KK %*% z) +
  if(!ch) {
    (1 - alpha) * t(z) %*% t(U) %*% M %*% U %*% z
  } else {
                 t(z) %*% t(U) %*% M %*% U %*% z
  }

  return(list(J = J,
              z = z,
              spline0 = spline0))
}
