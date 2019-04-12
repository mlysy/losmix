#' Genenerate random draws from the mNIX distribution.
#'
#' @param n Number of random draws (integer).
#' @param lambda Mean parameter.  A vector of length \code{p}, or a matrix of size \code{n x p}.
#' @param Omega Precision matrix parameter.  A matrix of size \code{p x p}, or an array of size \code{p x p x n}.
#' @param tau Scale parameter.  A scalar or a vector of length \code{n}.
#' @param nu Degrees-of-freedom parameter.  A scalar or a vector of length \code{n}.
#' @param y Optional response vector of length \code{N} for one subject.
#' @param X Optional covariate matrix of size \code{N x p} for one subject.
#' @return A list with elements \code{beta} and \code{sigma} of size \code{n x p} and length \code{n}, respectively.
#' @export
mnix_sim <- function(n, lambda, Omega, tau, nu, y, X) {
  # format inputs
  do_post <- !missing(y) || !missing(X)
  if(do_post) {
    # posterior distribution
    odata <- get_tmbdata(id = rep(1, length(y)), y = y, X = X)
    opars <- get_tmbpars(p = nrow(odata$Xtr),
                         lambda = lambda, Omega = Omega, tau = tau, nu = nu)
  } else {
    # arbitrary mnix distribution
    opars <- get_tmbpars(lambda = lambda, Omega = Omega, tau = tau, nu = nu)
    odata <- get_tmbdata(id = 1, y = 0, X = t(rep(0,nrow(opars$lambda))))
  }
  n <- c(opars$nOut, n)
  if(length(unique(n[n>1])) > 1) {
    stop("parameter size incompatible with n.")
  }
  opars$nOut <- max(n)
  odata <- c(list(model_name = "mNIX_sim", doPost = do_post), odata, opars)
  obj <- TMB::MakeADFun(data = odata, parameters = list(theta = 0),
                        silent = TRUE, DLL = "losmix_TMBExports")
  out <- obj$simulate()
  # reorder and resize outputs
  list(beta = t(out$beta), sigma = out$sigma)
}
