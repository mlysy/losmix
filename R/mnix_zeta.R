#' Calculate the mNIX normalizing constant.
#'
#' @template param-Omega
#' @template param-nu
#' @template param-tau
#' @return Vector of mNIX normalizing constants.
#' @details The mNIX normalzing constant (on the log scale) is defined as
#' \preformatted{
#' lgamma(nu/2) - .5 * ( nu * log(tau*nu/2) + log(|Omega|) )
#' }
#' @export
mnix_zeta <- function(Omega, nu, tau) {
  if(is.matrix(Omega)) Omega <- array(Omega, dim = c(dim(Omega),1))
  ldO <- apply(Omega, 3, function(om) solveV(om, ldV = TRUE)$ldV)
  lgamma(nu/2) - .5 * (nu * log(tau*nu/2) + ldO)
}
