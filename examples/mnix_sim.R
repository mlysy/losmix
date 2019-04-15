# helper functions

# generate random dataset with nSub subjects and p covariates
sim_data <- function(nSub, p) {
  N <- sample(2:10, nSub, replace = TRUE)
  M <- sum(N)
  X <- matrix(rnorm(M*p), M, p)
  y <- rnorm(M)
  id <- sample(rep(1:length(N), times = N))
  list(id = id, X = X, y = y)
}

# generate nPhi random hyperparameter sets on p covariates
sim_phi <- function(nPhi, p) {
  Phi <- replicate(n = nPhi, {
    list(lambda = rnorm(p), Omega = diag(rexp(p)),
         nu = rexp(1), tau = rexp(1))
  }, simplify = FALSE)
  # format form mnix_sim
  Phi <- unlist_bind(Phi,
                   name = c("lambda", "Omega", "nu", "tau"),
                   bind = c(rbind, cbind, c, c))
  Phi$Omega <- drop(array(Phi$Omega, dim = c(p, p, nPhi)))
  Phi
}

# sample from unrestricted mNIX distribution
n <- 5 # number of samples
p <- 2 # number of covariates
phi <- sim_phi(nPhi = 1, p = p) # generate hyperparameters
mnix_sim(n = n,
         lambda = phi$lambda, Omega = phi$Omega, nu = phi$nu, tau = phi$tau)

# vectorized calculations
Phi <- sim_phi(nPhi = n, p = p)
mnix_sim(n = n,
         lambda = Phi$lambda, Omega = Phi$Omega, nu = Phi$nu, tau = Phi$tau)


# sample from posterior mNIX distribution
# generate data
yX <- sim_data(nSub = n, p = p)

# single phi, single subject
mnix_sim(n = n,
         y = yX$y[yX$id == 1], X = yX$X[yX$id == 1,],
         lambda = phi$lambda, Omega = phi$Omega, nu = phi$nu, tau = phi$tau)

# single phi, multiple subjects
mnix_sim(n = n,
         y = yX$y, X = yX$X, id = yX$id,
         lambda = phi$lambda, Omega = phi$Omega, nu = phi$nu, tau = phi$tau)

# multiple phi, multiple subjects
mnix_sim(n = n,
         y = yX$y, X = yX$X, id = yX$id,
         lambda = Phi$lambda, Omega = Phi$Omega, nu = Phi$nu, tau = Phi$tau)
