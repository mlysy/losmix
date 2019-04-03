#--- basic helper functions ----------------------------------------------------

source("utils.R") # clean this up later

# log-determinant
ldet <- function(X) as.numeric(determinant(X)$modulus)

#--- simulation functions ------------------------------------------------------

# simulate data
sim_X <- function(N, p) matrix(rnorm(N*p), N, p)
sim_y <- function(N) rnorm(N)

# simulate random effects
sim_beta <- function(p) rnorm(p)
sim_sig2 <- function() rexp(1)
sim_theta <- function(p) list(beta = rnorm(p), sig2 = rexp(1))

# simulate hyperparameters
sim_Phi <- function(p) {
  lambda <- rnorm(p)
  Omega <- crossprod(matrix(rnorm(p^2), p, p))
  nu <- runif(1, 10, 20)
  tau <- rexp(1) + 5
  list(lambda = lambda, Omega = Omega, nu = nu, tau = tau)
}

#--- mniw methods --------------------------------------------------------------

# get sufficient statistics
get_suff <- function(y, X) {
  list(yy = crossprod(y)[1], Xy = crossprod(X, y),
       XX = crossprod(X), N = length(y))
}

# get hyperparameters of conjugate posterior
get_post <- function(suff, Phi) {
  list2env(suff, environment())
  list2env(Phi, environment())
  Ol <- Omega %*% lambda
  lOl <- crossprod(lambda, Ol)
  Ohat <- XX + Omega
  lhat <- c(solve(Ohat, Ol + Xy))
  nuhat <- N + nu
  that <- yy - crossprod(lhat, Ohat %*% lhat) + lOl
  that <- (that[1] + nu*tau)/nuhat
  list(lambda = lhat, Omega = Ohat, nu = nuhat, tau = that)
}

# normalizing constant
zeta <- function(Phi) {
  2 * lgamma(Phi$nu/2) - Phi$nu * log(Phi$tau*Phi$nu/2) - ldet(Phi$Omega)
}

#--- conversions between regular and computational basis -----------------------

# convert variance matrix to log-cholesky factor
var2lchol <- function(V) {
  C <- chol(V)
  diag(C) <- log(diag(C))
  C[upper.tri(C, diag = TRUE)]
}

# convert log-cholesky factor to variance matrix
lchol2var <- function(lC) {
  # determine size of matrix
  p <- (-1 + sqrt(1 + 8*length(lC)))/2
  C <- matrix(0, p,p)
  C[upper.tri(C, diag = TRUE)] <- lC
  diag(C) <- exp(diag(C))
  crossprod(C)
}

# convert hyperparameters to unconstained vector
Phi2vec <- function(Phi) {
  c(Phi$lambda, var2lchol(Phi$Omega), Phi$nu, Phi$tau)
}

# convert unconstrained vector to hyperparameter list
vec2Phi <- function(vphi) {
  # determine size of matrix
  p <- (-3 + sqrt(9 + 8*(length(vphi)-2)))/2
  lambda <- vphi[1:p]
  lC <- vphi[p+1:(p*(p+1)/2)]
  nu <- vphi[p*(p+3)/2+1]
  tau <- vphi[p*(p+3)/2+2]
  list(lambda = lambda, Omega = lchol2var(lC), nu = nu, tau = tau)
}

# convert hyperparameter list to TMB parameter list
Phi2par <- function(Phi) {
  list(lambda = Phi$lambda, Omega_LC = var2lchol(Phi$Omega),
       nu = Phi$nu, tau = Phi$tau)
}
