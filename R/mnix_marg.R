#' Marginal posterior of the mNIX hyperparameters as a \code{TMB::MakeADFun} object.
#'
#' @template param-id
#' @template param-y
#' @template param-X
#' @return A list as returned by \code{TMB::MakeADFun} corresponding to the \emph{negative} log-posterior of the marginal distribution of the mNIX hyperparameters (see \strong{Details}).
#'
#' @details In particular, the return list contains elements \code{fn}, \code{gr}, and \code{he} corresponding to the negative log-posterior and its first and second derivatives, computed analytically in C++ via automatic differentiation (AD).
#'
#' To evaluate these functions, the mNIX hyperparameters \code{Phi = (lambda, Omega, nu, tau)} must first be converted to an unconstrained vector, namely
#' \preformatted{
#' psi = (lambda, log_chol(Omega), log(nu), log(tau)),
#' }
#' where \code{link{log_chol}} is the log-Cholesky decomposition.
#'
#' The default prior is
#' \preformatted{
#' pi(Phi) \propto (nu * tau * |Omega|)^{-1},
#' }
#' which is equivalent to
#' \preformatted{
#' pi(psi) \propto prod( diag(chol(Omega))^(p-(1:p)) ).
#' }
#' @example examples/mnix_marg.R
#' @export
mnix_marg <- function(id, y, X) {
  # format inputs
  odata <- c(list(model_name = "mNIX_marg"), format_data(id = id, y = y, X = X))
  p <- nrow(odata$Xtr)
  opars <- list(lambda = rep(0, p), logC_Omega = log_chol(diag(p)),
                log_tau = 0, log_nu = 0)
  TMB::MakeADFun(data = odata, parameters = opars,
                 silent = TRUE, DLL = "losmix_TMBExports")
}
