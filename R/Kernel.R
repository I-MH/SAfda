#
#  this file is largely done and refactored; just need to check the docs (i and h), and write unit tests.
#
#' Kernel
#'
#' Given an input parameter `kerneltype`, returns the functional form of the kernel
#' as a function of inputs `i` and `h`: `x = i/h`.
#'
#' @param i fill in
#' @param h fill in
#' @param kerneltype Character input specifying type of kernel from list: `flat`, `Simple`, `Bartlett`, `flat_top` and `Parzen`
#'
#' @return Value of specified kernel for inputs `i` and `h`.
#' @export
#'
#' @examples Kernel(i = 1, h = 2, kerneltype = "Bartlett")
Kernel <- function(i, h, kerneltype) {
  stopifnot(kerneltype %in% c("flat", "Simple", "Bartlett", "flat_top", "Parzen"))
  stopifnot(is.numeric(i), is.numeric(h), is.finite(i / h))
  x <- i/h

  if (kerneltype == "flat") {
    return(1)
  } else if(kerneltype == "Simple") {
    return(0)
  } else if(kerneltype == "Bartlett") {
    return(1 - x)
  } else if (kerneltype == "flat_top") {
    if (x < 0.1) {
      return(1)
    } else {
      if (x >= 0.1 & x < 1.1) {
        return(1.1 - x)
      } else {
        return(0)
      }
    }
  } else if(kerneltype == "Parzen") {
    if (x < 1/2) {
      return(1 - 6 * x^2 + 6 * abs(x)^3)
    } else {
      return(2 * (1 - abs(x))^3)
    }
  }
}

