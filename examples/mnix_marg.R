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
# note that gradient information is included via the "attribute",
# to increase computational efficiency when function and gradient
# share calculations.
# this is the case for TMB functions, as calling nlp$gr() without
# arguments reuses calculations from the last call to nlp$fn.
ofun <- function(par) {
  out <- nlp$fn(par)
  # include gradient information via 'attribute'
  attr(out, "gradient") <- nlp$gr()
  out
}
opt <- nlm(p = nlp$par, # starting value
           f = ofun) # objective function

# gradient at the mode relative to its absolute value
nlp$gr(opt$estimate)/abs(opt$estimate)

## Approximate Bayesian Inference ##
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

# histograms of posterior distributions
# true hyperparameter values are vertical lines in red

# format data for plotting
iOmega <- cbind(rbind(1:p, 1:p), combn(1:p,2)) # unique elements of Omega
Phi_plot <- cbind(Phi_post$lambda, # lambda
                  apply(iOmega, 2,
                        function(ii) Phi_post$Omega[ii[1],ii[2],]), # Omega
                  Phi_post$nu, # nu
                  Phi_post$tau) # tau
# hyperparameter names
phi_names <- c(paste0("lambda[", 1:p, "]"),
               paste0("Omega[", iOmega[1,], iOmega[2,], "]"),
               "nu", "tau")
# true values
phi_true <- c(lambda, Omega[t(iOmega)], nu, tau)

# create plot
par(mfrow = c(3,4), mar = c(2,2,4,.5))
for(ii in 1:ncol(Phi_plot)) {
  # approximate posterior
  hist(Phi_plot[,ii], breaks = 40, xlab = "", ylab = "",
       main = parse(text = paste0("hat(p)(", phi_names[ii],
                                  "*\" | \"*bold(Y),bold(X))")))
  # true parameter value
  abline(v = phi_true[ii], col = "red", lwd = 2)
}
# legend
plot(0, type = "n", xlim = c(0,1), ylim = c(0,1),
     xlab = "", ylab = "", axes = FALSE)
legend("bottom", inset = .05,
       legend = c("Posterior Distribution", "True Hyperparameter Value"),
       lwd = c(NA, 2), pch = c(22, NA), seg.len = 1.5,
       col = c("black", "red"), bg = c("white", NA), cex = .85)

## Inference for Random Effects ##

iSub <- sample(nSub, 1) # pick a subject at random
# data for subject i
Xi <- X[id == iSub,]
yi <- y[id == iSub]
# true parameter values for subject i
thetai_true <- c(Theta$beta[iSub,], sigma = Theta$sigma[iSub])

# method 1: using mnix_sim
set.seed(1)
Thetai_post <- mnix_sim(n = npost, y = yi, X = Xi,
                        lambda = Phi_post$lambda, Omega = Phi_post$Omega,
                        nu = Phi_post$nu, tau = Phi_post$tau)

# method 2: using mnix_marg
nlpi <- mnix_marg(y = yi, X = Xi) # posterior distribution given just subject i
set.seed(1)
Thetai_post2 <- apply(Psi_post, 1, function(psi) nlpi$simulate(psi))
# reshape data
Thetai_post2 <- unlist_bind(Thetai_post2,
                            name = c("beta", "sigma"),
                            bind = c(cbind, c))
Thetai_post2$beta <- t(Thetai_post2$beta) # beta needs to be transposed

# check that they are the same
all.equal(Thetai_post, Thetai_post2)

# format data for plotting
Thetai_plot <- cbind(Thetai_post$beta, # beta
                     Thetai_post$sigma) # sigma
# parameter names
thetai_names <- c(paste0("beta[i", 1:p, "]"), "sigma[i]")
# true parameter values for subject i
thetai_true <- c(Theta$beta[iSub,], sigma = Theta$sigma[iSub])

# create plot
par(mfrow = c(2,2), mar = c(2,2,4,.5))
for(ii in 1:ncol(Thetai_plot)) {
  # approximate posterior
  hist(Thetai_plot[,ii], breaks = 40, xlab = "", ylab = "",
       main = parse(text = paste0("hat(p)(", thetai_names[ii],
                                  "*\" | \"*bold(Y),bold(X))")))
  # true parameter value
  abline(v = thetai_true[ii], col = "red", lwd = 2)
}
