#--- basic helper functions ----------------------------------------------------

source("utils.R") # clean this up later

# log-determinant
ldet <- function(X) as.numeric(determinant(X)$modulus)

# sort list on element names
sort_names <- function(list) {
  if(!is.null(nm <- names(list))) list[sort(nm)] else list
}

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

# simulate ids
sim_id <- function(N) {
  M <- length(N)
  lab <- sample(1000, M)
  id <- data.frame(id = sample(rep(lab, times = N)),
                   order = NA)
  for(jj in 1:M) id$order[id$id==lab[jj]] <- cumsum(c(0,N))[jj]+1:N[jj]
  id$id <- factor(id$id, levels = lab)
  id
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

#--- conversions between regular and computational basis -----------------------


# convert hyperparameter list to TMB parameter list
Phi2par <- function(Phi) {
  list(lambda = Phi$lambda, logC_Omega = losmix::log_chol(Phi$Omega),
       log_nu = log(Phi$nu), log_tau = log(Phi$tau))
}
