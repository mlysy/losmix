#' Calculate the parameters of the mNIX conjugate posterior distribution.
#'
#' @template param-id
#' @template param-y
#' @template param-X
#' @param lambda Prior mean parameter.  A vector of length \code{p}, or a matrix of size \code{n x p}.
#' @param Omega Prior precision matrix.  A matrix of size \code{p x p}, or an array of size \code{p x p x n}.
#' @param nu Prior degrees-of-freedom.  A scalar or a vector of length \code{n}.
#' @param tau Prior scale parameter.  A scalar or a vector of length \code{n}.
#' @return A list with elements \code{lambda}, \code{Omega}, \code{nu}, and \code{tau} of posterior parameters, consisting of a matrix of size \code{n x p}, and array of size \code{p x p x n}, and two vectors of length \code{n}, respectively.
#' @export
mnix_post <- function(id, y, X, lambda, Omega, nu, tau) {
  # format inputs
  odata <- format_data(id = id, y = y, X = X)
  p <- nrow(odata$Xtr)
  opars <- format_pars(p = p,
                       lambda = lambda, Omega = Omega, tau = tau, nu = nu)
  n <- c(opars$nOut, length(odata$nObs))
  if(length(unique(n[n>1])) > 1) {
    stop("parameter size incompatible with length(unique(id)).")
  }
  opars$nOut <- max(n)
  odata <- c(list(model_name = "mNIX_post"), odata, opars)
  obj <- TMB::MakeADFun(data = odata, parameters = list(theta = 0),
                        silent = TRUE, DLL = "losmix_TMBExports")
  out <- obj$simulate()
  # relabel and reorder outputs
  out_names <- c(lambda = "lambda_hat", Omega = "Omega_hat",
                 tau = "tau_hat", nu = "nu_hat")
  out <- out[out_names]
  names(out) <- names(out_names)
  out$lambda <- t(out$lambda)
  out$Omega <- array(out$Omega, dim = c(p,p,opars$nOut))
  out
}
