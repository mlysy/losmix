# simulate data

p <- 3 # number of covariates
nsub <- 10 # number of subjects
p <- 3 # number of covariates
N <- sample(20:50, nsub, replace = TRUE) # number of

# hyperparameters
lambda <- rnorm(p)
Omega <- diag(p)
nu <- 5 * runif(1) + 3
tau <- rexp(1)

# parameters
Theta <- mnix_sim(nsub, lambda = lambda, Omega = Omega, nu = nu, tau = tau)
