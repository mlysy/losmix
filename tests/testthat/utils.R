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

# density of mnix distribution.
# x is a vector.
# v is a scalar.
# Phi is a list with elements lambda, Omega, nu, tau.
dmnix <- function(x, v, Phi, log = FALSE) {
  ans <- losmix::dsichisq(v, nu = Phi$nu, tau = Phi$tau, log = TRUE)
  ans <- ans + losmix::dmvn(x, mu = Phi$lambda, Sigma = v * solve(Phi$Omega),
                            log = TRUE)
  if(!log) ans <- exp(ans)
  ans
}

# random draw from the mnix distribution
rmnix <- function(lambda, Omega, nu, tau) {
  sigma <- sqrt(losmix::rsichisq(1, nu, tau))
  z <- rnorm(length(lambda), sd = sigma)
  beta <- backsolve(chol(Omega), z) + lambda
  list(beta = c(beta), sigma = sigma)
}
