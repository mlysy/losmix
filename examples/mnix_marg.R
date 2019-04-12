# simulate data

set.seed(1234) # reproducible results

p <- 3 # number of covariates
nSub <- 50 # number of subjects
N <- sample(20:50, nSub, replace = TRUE) # observations per subject

# hyperparameters
lambda <- rnorm(p)
Omega <- diag(p)
nu <- 5 * runif(1) + 3
tau <- rexp(1)

# parameters
Theta <- mnix_sim(nSub, lambda = lambda, Omega = Omega, nu = nu, tau = tau)

# data
# design matrices
X <- lapply(N, function(n) matrix(rnorm(n*p), n, p))
# response vectors
y <- lapply(1:nSub, function(ii) {
  rnorm(N[ii], mean = c(X[[ii]] %*% Theta$beta[ii,]), sd = Theta$sigma[ii])
})
# convert lists to single matrix/vector
X <- do.call(rbind, X)
y <- do.call(c, y)
# subject identifiers
id <- rep(1:nSub, times = N)

# marginal posterior distribution
nlp <- mnix_marg(id = id, y = y, X = X)

# mode-finding using quasi-Newton algorithm in nlm

# objective function
obj_fun <- function(par) {
  out <- nlp$fn(par)
  # include gradient information via attribute
  attr(out, "gradient") <- nlp$gr(par)
  out
}
opt <- nlm(p = nlp$par, # starting value
           f = obj_fun)

# gradient at the mode relative to its absolute value
nlp$gr(opt$estimate)/abs(opt$estimate)

# ## Approximate Bayesian Inference ##
#
# On transformed scale:
#
# - posterior mean is mode of nlp
# - posterior variance is inverse of hessian of nlp
#
# On original scale:
#
# - simulate from multivariate normal on the transformed scale
# - convert output to original scale

# transformed scale simulation
npost <- 1e4 # number of posterior draws
psi_mean <- opt$estimate # posterior mean
psi_var <- solveV(nlp$he(opt$estimate))
Psi_post <- rmvn(n = npost, mu = psi_mean, Sigma = psi_var)

# convert to original scale
Phi_post <- ivec_phi(Psi_post)

# plot (approximate) posterior distributions
# true parameter values are plotted in red
par(mfrow = c(2,4))
# lambda
for(ii in 1:p) {
  hist(Phi_post$lambda[,ii],
       freq = FALSE, breaks = 50, xlab = "",
       main = parse(text = paste0("lambda[", ii, "]")))
  abline(v = lambda[ii], col = "red")
}
# Omega
for(ii in list(c(1,1), c(1,2), c(2,2))) {
  hist(Phi_post$Omega[ii[[1]],ii[[2]],],
       freq = FALSE, breaks = 50, xlab = "",
       main = parse(text = paste0("Omega[", ii[[1]], ii[[2]], "]")))
  abline(v = Omega[ii[[1]],ii[[2]]], col = "red")
}
# nu
hist(Phi_post$nu,
     freq = FALSE, breaks = 50, xlab = "",
     main = parse(text = "nu"))
abline(v = nu, col = "red")
# tau
hist(Phi_post$tau,
     freq = FALSE, breaks = 50, xlab = "",
     main = parse(text = "tau"))
abline(v = tau, col = "red")
