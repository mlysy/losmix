#' Scale-invariant variance matrix prior on the log-Cholesky scale.
#'
#' @template param-Omega
#' @return A log-prior vector of length \code{n}.
#' @details For a given variance (or precision) matrix \code{Omega}, the scale-invariant prior is
#' \preformatted{
#' pi(Omega) \propto 1/|Omega|
#' }
#' in the sense that the change of variables \code{Psi = A Omega t(A)} results in the prior \code{pi(Psi) = 1/|Psi|}.
#'
#' On the log-Cholesky scale defined by \code{\link{log_chol}}, the scale-invariant prior is given by
#' \preformatted{
#' pi(logC) \propto prod( diag(chol(Omega))^(p-(1:p)) )
#' }
#' @export
lchol_prior <- function(Omega) {
  if(is.matrix(Omega)) Omega <- array(Omega, dim = c(dim(Omega),1))
  omega <- apply(Omega, 3, function(om) log(diag(chol(om))))
  omega <- as.matrix(omega)
  p <- nrow(omega)
  colSums((p-(1:p)) * omega)
}
