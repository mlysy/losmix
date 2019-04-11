#' Compute the log-Cholesky decomposition and its inverse.
#'
#' @param V A variance matrix of size \code{p x p}.
#' @param logC The log-Cholesky decomposition, represented as a vector of length \code{p*(p+1)/2} (see \strong{Details}).
#' @return For \code{log_chol}, the log-Cholesky decomposition (a vector of length \code{p*(p+1)/2}).  For \code{ilog_chol}, a variance matrix of size \code{p x p}.
#' @details The log-Cholesky decomposition of a variance matrix \code{V} is the upper triangular matrix, of which the non-zero elements are concatenated in column-major order.  Namely, the calculation is given by
#' \preformatted{
#' logC <- chol(V)
#' diag(logC) <- log(diag(logC))
#' logC <- logC[upper.tri(logC,diag=TRUE)]
#' }
#' @name log_chol
#' @export
log_chol <- function(V) {
  logC <- chol(V)
  diag(logC) <- log(diag(logC))
  logC[upper.tri(logC, diag = TRUE)]
}

#' @rdname log_chol
#' @export
ilog_chol <- function(logC) {
  # determine size of matrix
  p <- (-1 + sqrt(1 + 8*length(logC)))/2
  C <- matrix(0, p,p)
  C[upper.tri(C, diag = TRUE)] <- logC
  diag(C) <- exp(diag(C))
  crossprod(C)
}
