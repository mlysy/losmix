#' Sufficient statistics for the mNIX distribution.
#'
#' @template param-id
#' @template param-y
#' @template param-X
#' @export
mnix_suff <- function(id, y, X) {
  # convert to TMB input
  odata <- c(list(model_name = "mNIX_suff"), get_tmbdata(id = id, y = y, X = X))
  opars <- list(theta = 0)
  obj <- TMB::MakeADFun(data = odata, parameters = opars,
                        silent = TRUE, DLL = "losmix_TMBExports")
  out <- obj$simulate()
  out[c("yy", "Xy", "XX")]
}

