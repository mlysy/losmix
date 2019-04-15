#' Convert between hyperparameter list and vector representation.
#'
#' @param Phi Hyperparameters as a list with elements \code{lambda}, \code{Omega}, \code{nu}, and \code{tau}.
#' @param psi Hyperparameters as an unconstrained numeric vector, or a matrix where each row is one hyperparameter set (see \strong{Details}).
#' @name vec_phi
#' @return For \code{vec_phi}, an unconstrained vector or matrix.  For \code{ivec_phi}, a list with elements \code{lambda}, \code{Omega}, \code{nu}, and \code{tau}.
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
  if(!is.matrix(psi)) psi <- t(psi)
  # determine size of matrix
  p <- (-3 + sqrt(9 + 8*(ncol(psi)-2)))/2
  lambda <- psi[,1:p]
  logC <- psi[,p+1:(p*(p+1)/2),drop=FALSE]
  ilc <- TMB::MakeADFun(data = list(model_name = "ilog_chol",
                                    logC = t(logC)),
                        parameters = list(theta = 0),
                        silent = TRUE, DLL = "losmix_TMBExports")
  Omega <- ilc$simulate()$V
  ## Omega2 <- apply(logC, 1, ilog_chol)
  Omega <- drop(array(Omega, dim = c(p, p, nrow(psi))))
  nu <- exp(psi[,p*(p+3)/2+1])
  tau <- exp(psi[,p*(p+3)/2+2])
  list(lambda = lambda, Omega = Omega, nu = nu, tau = tau)
}
