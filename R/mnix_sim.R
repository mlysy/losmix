#' Genenerate random draws from the mNIX distribution.
#'
#' @param n Number of random draws (integer).
#' @template param-lambda
#' @template param-Omega
#' @template param-nu
#' @template param-tau
#' @param y Optional response vector of length \code{N} for one subject.
#' @param X Optional covariate matrix of size \code{N x p} for one subject.
#' @return A list with elements \code{beta} and \code{sigma} of size \code{n x p} and length \code{n}, respectively.
#' @export
mnix_sim <- function(n, lambda, Omega, nu, tau, y, X) {
  # format inputs
  do_post <- !missing(y) || !missing(X)
  if(do_post) {
    # posterior distribution
    odata <- format_data(y = y, X = X)
    opars <- format_pars(p = nrow(odata$Xtr),
                         lambda = lambda, Omega = Omega, tau = tau, nu = nu)
  } else {
    # arbitrary mnix distribution
    opars <- format_pars(lambda = lambda, Omega = Omega, tau = tau, nu = nu)
    odata <- format_data(y = 0, X = t(rep(0,nrow(opars$lambda))))
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
