#' Convert between hyperparameter list and vector representation.
#'
#' @param Phi Hyperparameters as a list with elements \code{lambda}, \code{Omega}, \code{nu}, and \code{tau}.
#' @param psi Hyperparameters as an unconstrained numeric vector (see \strong{Details}).
#' @name vec_phi
#' @return For \code{vec_phi}, an unconstrained vector.  For \code{ivec_phi}, a list with elements \code{lambda}, \code{Omega}, \code{nu}, and \code{tau}.
#' @details The unconstrained hyperparameter representation is given by
#' \preformatted{
#' psi <- c(lambda, log_chol(Omega), log(nu), log(tau))
#' }
#' where \code{\link{log_chol}} is the log-Cholesky decomposition.
#' @export
# convert hyperparameters to unconstained vector
vec_phi <- function(Phi) {
  c(Phi$lambda, log_chol(Phi$Omega), log(Phi$nu), log(Phi$tau))
}

#' @rdname vec_phi
#' @export
ivec_phi <- function(psi) {
  # determine size of matrix
  p <- (-3 + sqrt(9 + 8*(length(psi)-2)))/2
  lambda <- psi[1:p]
  logC <- psi[p+1:(p*(p+1)/2)]
  nu <- exp(psi[p*(p+3)/2+1])
  tau <- exp(psi[p*(p+3)/2+2])
  list(lambda = lambda, Omega = ilog_chol(logC), nu = nu, tau = tau)
}
