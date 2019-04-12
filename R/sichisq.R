#' Density and random sampling from the scaled inverse-chi-square distribution.
#' @name sichisq
#'
#' @param x Observations (positive scalar or vector).
#' @param n Integer number of random draws to generate.
#' @param nu Shape parameter (positive scalar or vector).
#' @param tau Scale parameter (positive scalar or vector).
#' @param log Logical; whether to calculate the density on the log scale.
#' @details Arguments \code{x}, \code{nu}, and \code{tau} are recycled as needed.
#' @export
rsichisq <- function(n, nu, tau) {
  1/rgamma(n, shape = nu/2, rate = nu*tau/2)
}

#' @rdname sichisq
#' @export
dsichisq <- function(x, nu, tau, log = FALSE) {
  ans <- dgamma(1/x, shape = nu/2, rate = nu*tau/2, log = TRUE) - 2*log(x)
  if(!log) ans <- exp(ans)
  ans
}
