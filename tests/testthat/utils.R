#--- utility functions ---------------------------------------------------------

# Density of the Inverse-Gamma distribution.
#
# @param x Vector of quantiles at which to evaluate the density.
# @param shape Vector of shape parameters.
# @param scale Vector of scale parameters.
# @param log Logical; whether or not to return the density of the log scale.
# @return Vector of density evaluations.
# @details The Inverse-Gamma distribution here with shape parameter \code{a} and scale parameter \code{b} has PDF
# \preformatted{
# f(x | a, b) = cst * x^(a + 1) * exp(-b/x),   cst = b^a/Gamma(a).
# }
# The distribution is proper only when \code{a,b > 0}, but these values are unrestricted here since improper priors can be useful for Bayesian inference.  When either \code{a,b <= 0} the function sets \code{cst = 1} in the expression above.
dinvgamma <- function(x, shape, scale, log = FALSE) {
  if(all(shape > 0 & scale > 0)) {
    cst <- shape * log(scale) - lgamma(shape)
  } else cst <- 0
  ld <- cst - (shape + 1) * log(x) - scale/x
  if(!log) ld <- exp(ld)
  ld
}

# Random draws from Inverse-Gamma distribution.
rinvgamma <- function(n, shape, scale) {
  1/rgamma(n, shape = shape, rate = scale)
}

# Random draws from a multivariate normal with mean mu and variance V.
rmvn <- function(n, mu, V) {
  p <- nrow(V)
  x <- matrix(rnorm(p*n), n, p) %*% chol(V)
  t(t(x) + mu)
}

# Calcuate the density of a multivariate normal with mean mu and variance V.
dmvn <- function(x, mu, V, log = FALSE) {
  z <- x-mu
  IP <- solveV(V, z, ldV = TRUE)
  ld <- -.5 * (sum(z*IP$y) + IP$ldV + length(mu) * log(2*pi))
  if(!log) ld <- exp(ld)
  ld
}

# Solve method for variance matrices.
#
# @param V Variance matrix
# @param x Optional vector or matrix for which to solve system of equations.  If missing calculates inverse matrix.
# @param ldV Optionally compute log determinant as well.
# @return Matrix solving system of equations and optionally the log-determinant.
# @details This function is faster and more stable than \code{solve} when \code{V} is known to be positive-definite.
solveV <- function(V, x, ldV = FALSE) {
  C <- chol(V)
  if(missing(x)) x <- diag(nrow(V))
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }
  ans
}

# density of scaled-inverse chi-square
# for testing, always use form which is least likely to contain errors
# in this case, this is change-of-variables y = 1/x with gamma distribution.
dsichisq <- function(x, nu, tau, log = FALSE) {
  ans <- dgamma(1/x, shape = nu/2, rate = nu*tau/2, log = TRUE) - 2*log(x)
  if(!log) ans <- exp(ans)
  ans
}

# density of mnix distribution.
# x is a vector.
# v is a scalar.
# Phi is a list with elements lambda, Omega, nu, tau.
dmnix <- function(x, v, Phi, log = FALSE) {
  ans <- dsichisq(v, nu = Phi$nu, tau = Phi$tau, log = TRUE)
  ans <- ans + dmvn(x, mu = Phi$lambda, V = v * solve(Phi$Omega),
                    log = TRUE)
  if(!log) ans <- exp(ans)
  ans
}
