#' Solve method for variance matrices.
#'
#' @param V Symmetric positive-definite (i.e., variance) matrix.
#' @param b Optional vector or matrix for which to solve the system of equations \code{V x = b}.  If missing calculates the inverse of \code{V}.
#' @param ldV Optionally compute the log determinant as well.
#' @return The solution of the system of equations, or a list with elements \code{x} and \code{ldV} returning the log-determinant as well.
#' @details This function is faster and more stable than \code{base::solve} when \code{V} is known to be positive-definite.
#' @export
solveV <- function(V, b, ldV = FALSE) {
  C <- chol(V)
  if(missing(b)) b <- diag(nrow(V))
  ans <- backsolve(r = C, x = backsolve(r = C, x = b, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(x = ans, ldV = ldV)
  }
  ans
}


#' Format data for \pkg{TMB} calculations.
#'
#' @template param-y
#' @template param-X
#' @template param-id
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{y}}{The response vector \code{y}, reordered such that all observations for a given subject are consecutive, and converted to an \code{n x 1} column matrix.}
#'   \item{\code{Xtr}}{The transpose of the covariate matrix, reordered the same way as \code{y}, having size \code{p x n}.}
#'   \item{\code{iStart}}{An integer vector of length \code{nsub}, specifying the starting index of the observations for each subject, using C++ indexing (i.e., starting at zero).}
#'   \item{\code{nObs}}{An integer vector of length \code{nsub}, specifying the number of observations per subject.}
#' }
#' @export
format_data <- function(y, X, id) {
  ntot <- length(y) # total number of observations
  if(missing(id)) id <- rep(1, ntot)
  if(!is.matrix(X)) X <- as.matrix(X)
  if(length(id) != ntot || nrow(X) != ntot) {
    stop("id, y, and X have incompatible dimensions.")
  }
  isub <- tapply(1:ntot, id, c, simplify = FALSE) # indices per subject
  # convert id to (istart, Ni)
  Ni <- sapply(isub, length)
  istart <- cumsum(c(0, Ni[-length(Ni)]))
  # reorder y and X
  y <- do.call(c, lapply(isub, function(ind) y[ind]))
  names(y) <- NULL
  X <- do.call(rbind, lapply(isub, function(ind) X[ind,,drop=FALSE]))
  list(y = as.matrix(y), Xtr = t(X), iStart = istart, nObs = Ni)
}

# format tmb parameters.  throw errors if there are problems.
format_pars <- function(p, nData, nSim, lambda, Omega, tau, nu) {
  if(is_narray(lambda, 1:2)) {
    lambda <- if(is.matrix(lambda)) t(lambda) else matrix(lambda)
    if(missing(p)) p <- nrow(lambda)
    if(nrow(lambda) != p) stop("dimensions of lambda are incompatible with p.")
  } else {
    stop("lambda must be a numeric vector or matrix.")
  }
  if(is_narray(Omega, 2:3)) {
    if(!all(dim(Omega)[1:2] == p)) stop("dimensions of Omega are incompatible with p.")
    Omega <- matrix(Omega, nrow = p)
  } else {
    stop("Omega must be a numeric matrix or array.")
  }
  if(!is_narray(nu, 1)) {
    stop("nu must be a numeric scalar or matrix.")
  }
  if(!is_narray(tau, 1)) {
    stop("tau must be a numeric scalar or matrix.")
  }
  # determine n
  n <- c(ncol(lambda), ncol(Omega)/p, length(nu), length(tau))
  if(length(unique(n[n>1])) > 1) {
    stop("lambda, Omega, nu, and tau have incompatible lengths.")
  }
  list(nOut = max(n), lambda = lambda, Omega = Omega, nu = nu, tau = tau)
}

# check if x is a numeric array with length(dim(x)) %in% dims
is_narray <- function(x, dims) {
  is.numeric(x) && (length(dim(as.array(x))) %in% dims)
}
