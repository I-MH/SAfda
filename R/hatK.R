#
#  this file is largely done and refactored; just need to check the docs, and write unit tests,
#  and put in sanity checks at the start
#
#' hatK
#'
#' Computes the auxiliary K, as the ratio of each eigenvalue to its predecessor;
#' returns the index and value of the minimum of such quotients.
#'
#' @param evalues  Numeric vector of eigenvalues
#'
#' @return list
#' hat_k  Index of the best k
#' ratio_k  The value of the best k
#' @export
#'
#' @examples rnorm(1)
hatK <- function(evalues){
  k_aux <- evalues[-1] / evalues[-length(evalues)]
  index_k <- which.min(k_aux)
  return(list(hat_k = index_k, ratio = k_aux))
}
