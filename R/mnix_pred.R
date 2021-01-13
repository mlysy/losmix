#' Genenerate random draws from the mNIX predictive distribution.
#'
#' @param n Number of random draws (integer).
#' @param Xpred Covariate matrix of size \code{n x p} to use for the predictive distribution.
#' @template param-lambda
#' @template param-Omega
#' @template param-nu
#' @template param-tau
#' @param y Optional response vector of length \code{N} (see \strong{Details}).
#' @param X Optional covariate matrix of size \code{N x p} (see \strong{Details}).
#' @param id Optional subject identifier vector of length \code{N} with \code{nsub <= N} distinct elements (see \strong{Details}).
#' @return A vector of \code{n} draws from the predictive distribution.
#'
#' @details If \code{y}, \code{X} and \code{id} are unspecified, then \code{(beta, sigma)} are drawn from \code{mNIX(lambda, Omega, nu, tau)}.  Otherwise, they are drawn from the posterior mNIX distribution \code{p(beta, sigma | y, X, id, lambda, Omega, nu, tau)}.  In this case, the number of subjects \code{nsub} must be \code{n} or \code{1}.
#' @export
mnix_pred <- function(n, Xpred, lambda, Omega, nu, tau, y, X, id) {
  # not "maximally" effficient, but will do.
  # so first sample from random effects posterior
  Theta <-  mnix_sim(n, lambda = lambda, Omega = Omega, nu = nu, tau = tau,
                     y = y, X = X, id = id)
  # now sample from predictive distribution
  rnorm(n, mean = rowSums(Xpred * Theta$beta), sd = Theta$sigma)
}
