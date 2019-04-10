# convert (id, y, X) to TMB input:
# 1. group y and X by consecutive ids
# 2. convert id to istart: starting point, and Ni: number of observations per subject.
get_tmbdata <- function(id, y, X) {
  if(!is.matrix(X)) X <- as.matrix(X)
  ntot <- length(id) # total number of observations
  if(length(y) != ntot || nrow(X) != ntot) {
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
  list(y = y, X = X, iStart = istart, nObs = Ni)
}

# format tmb parameters.  throw errors if there are problems.
get_tmbpars <- function(p, nData, lambda, Omega, tau, nu) {
  if(is_narray(lambda, 1:2)) {
    lambda <- if(is.matrix(lambda)) t(lambda) else matrix(lambda)
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
  n <- c(max(n), nData)
  if(length(unique(n[n>1])) > 1) {
    stop("parameter size incompatible with length(unique(id)).")
  }
  list(nPost = max(n), lambda = lambda, Omega = Omega, nu = nu, tau = tau)
}

# check if x is a numeric array with length(dim(x)) %in% dims
is_narray <- function(x, dims) {
  is.numeric(x) && (length(dim(as.array(x))) %in% dims)
}
