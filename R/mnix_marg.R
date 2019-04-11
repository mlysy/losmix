#' Marginal posterior of the mNIX hyperparameters as a \code{TMB::MakeADFun} object.
#'
#' @template param-id
#' @template param-y
#' @template param-X
#' @export
mnix_marg <- function(id, y, X) {
  # format inputs
  odata <- c(list(model_name = "mNIX_marg"), get_tmbdata(id = id, y = y, X = X))
  p <- ncol(odata$X)
  opars <- list(lambda = rep(0, p), logC_Omega = log_chol(diag(p)),
                log_tau = 0, log_nu = 0)
  TMB::MakeADFun(data = odata, parameters = opars,
                 silent = TRUE, DLL = "losmix_TMBExports")
}
