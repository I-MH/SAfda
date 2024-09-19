#' hatK
#'
#' Computes the auxiliary K, as the ratio of each eigenvalue to its predecessor;
#' returns the index and the ratio values.
#'
#' @param evalues  Numeric vector of eigenvalues
#'
#' @return 
#'  A list with: 
#' \itemize{
#' \item \code{hat_k}:  Index of the best k
#' \item \code{ratio_k}:  values of the ratio
#' }

#' @export
#'
#' @examples 
#' x <- rnorm(10)
#' hatK(x)
hatK <- function(evalues){
  k_aux <- evalues[-1] / evalues[-length(evalues)]
  index_k <- which.min(k_aux)
  return(list(hat_k = index_k, ratio = k_aux))
}
