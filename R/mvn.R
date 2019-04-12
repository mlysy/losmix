#' Density and random sampling for the multivariate normal distribution.
#' @name mvn
#'
#' @param x Mutivariate observations.  A vector or length \code{n} or a matrix of size \code{n x p}.
#' @param n Integer number of random draws.
#' @param mu Mean vector of length \code{p}.
#' @param Sigma Variance matrix of size \code{p x p}.
#' @param log Logical; whether to calculate the density on the log scale.
#' @return For \code{rmvn}, an \code{n x p} matrix of random draws.  For \code{dmvn} a vector of \code{n} density evaluations.
#' @export
rmvn <- function(n, mu, Sigma) {
  p <- nrow(Sigma)
  x <- matrix(rnorm(p*n), n, p) %*% chol(Sigma)
  t(t(x) + mu)
}

#' @rdname mvn
#' @export
dmvn <- function(x, mu, Sigma, log = FALSE) {
  if(!is.matrix(x)) x <- t(x)
  z <- t(x)-mu
  IP <- solveV(Sigma, z, ldV = TRUE)
  ld <- -.5 * (colSums(z*IP$x) + IP$ldV + length(mu) * log(2*pi))
  if(!log) ld <- exp(ld)
  ld
}
