#' Sufficient statistics for the mNIX distribution.
#'
#' @template param-id
#' @template param-y
#' @template param-X
#' @return A list with elements:
#' \describe{
#'   \item{\code{yy}}{A vector of length \code{nsub} consisting of \code{t(y) y} for each subject.}
#'   \item{\code{Xy}}{A matrix of size \code{p x nsub} consisting of \code{t(X) y} for each subject.}
#'   \item{\code{XX}}{A matrix of size \code{p x nsub^2} consisting of \code{t(X) X} for each subject.}
#' }
#' @export
mnix_suff <- function(id, y, X) {
  # convert to TMB input
  odata <- c(list(model_name = "mNIX_suff"), format_data(id = id, y = y, X = X))
  opars <- list(theta = 0)
  obj <- TMB::MakeADFun(data = odata, parameters = opars,
                        silent = TRUE, DLL = "losmix_TMBExports")
  out <- obj$simulate()
  out[c("yy", "Xy", "XX")]
}

