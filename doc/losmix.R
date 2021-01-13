## ----setup, include = FALSE---------------------------------------------------
suppressMessages({
  require(losmix)
  require(TMB)
})
source("format.R")
knitr::opts_chunk$set(comment = NA)

## ----ex1_sim------------------------------------------------------------------
require(losmix)

set.seed(1234) # reproducible results

p <- 3 # number of covariates
nSub <- 50 # number of subjects
N <- sample(20:50, nSub, replace = TRUE) # number of observations per subject

# hyperparameters
lambda <- rnorm(p)
Omega <- diag(p)
nu <- 5 * runif(1) + 3
tau <- rexp(1)

# parameters: generate them from an mNIX distribution
Theta <- mnix_sim(nSub, lambda = lambda, Omega = Omega, nu = nu, tau = tau)

# data: generate vector/matrix for each subject,
# then merge into single vector/matrix
# covariate matrices
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

## ----ex1_nlp------------------------------------------------------------------
# marginal posterior distribution
nlp <- mnix_marg(id = id, y = y, X = X)
names(nlp)

## ----ex1_fit, warning = FALSE-------------------------------------------------
# objective function
ofun <- function(par) {
  out <- nlp$fn(par)
  # include gradient information via 'attribute'
  attr(out, "gradient") <- nlp$gr()
  out
}

# optimization
opt <- nlm(p = nlp$par, # starting value (losmix picks a reasonable default)
           f = ofun) # objective function
opt$code # code == 1 means that 'nlm' thinks it converged

## ----ex1_grad_check-----------------------------------------------------------
disp <- rbind(est = opt$estimate, # potential solution
              grad = opt$gradient, # gradient at the potential solution
              rel = opt$gradient/abs(opt$estimate)) # relative size
signif(disp, 2)

## ----ex1_hess-----------------------------------------------------------------
psi_mean <- opt$estimate # (approximate) posterior mean of p(psi | Y, X)
psi_var <- solveV(nlp$he(opt$estimate)) # (approximate) posterior variance

## ----ex1_mpost, fig.width = 7, fig.height = 5---------------------------------
# Step 2: sample from p_hat(psi | Y, X)

npost <- 1e4 # number of posterior draws
Psi_post <- rmvn(n = npost, mu = psi_mean, Sigma = psi_var)

# Step 3: convert to sample from p_hat(phi | Y, X)
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

## ----ex1_rxpost, fig.width = 7, fig.height = 5--------------------------------
# inference for random effects

iSub <- sample(nSub, 1) # pick a subject at random

# data for subject i
Xi <- X[id == iSub,]
yi <- y[id == iSub]

# sample from p(thetai | y, X)
Thetai_post <- mnix_sim(npost,
                        lambda = Phi_post$lambda, Omega = Phi_post$Omega,
                        nu = Phi_post$nu, tau = Phi_post$tau,
                        y = yi, X = Xi)

# plot (approximate) posterior distributions
# true parameter values are plotted in red

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

## ----ex2_dir, echo = 2, results = "hide"--------------------------------------
mxt_file <- system.file("include", "losmix", "ModelExt.cpp", package = "losmix")
system.file("include", "losmix", "ModelExt.cpp", package = "losmix")

## ---- echo = FALSE, results = "asis"------------------------------------------
cat("```cpp", readLines(mxt_file), "```", sep = "\n")

## ----ex2_compile, eval = FALSE------------------------------------------------
#  require(TMB)
#  
#  model_name <- "ModelExt"
#  # instruct TMB where to find losmix library
#  include_path <- system.file("include", package = "losmix")
#  # compile and load model
#  TMB::compile(paste0(model_name, ".cpp"),
#               PKG_CXXFLAGS = paste0('-I"', include_path, '"'))
#  dyn.load(TMB::dynlib(model_name))

## ----ex2_data, fig.width = 7, fig.height = 4----------------------------------
nSub <- 10 # number of subjects
N <- sample(20:50, nSub, replace = TRUE) # observations per subject

# hyperparameters
gamma <- runif(1) * .2
lambda <- rnorm(2)
Omega <- crossprod(matrix(rnorm(4), 2, 2))/10
nu <- runif(1, 1, 2)*50
tau <- rexp(1)/100

# covariates
X <- lapply(N, function(N) {
  cbind(t = 1:N, x = runif(1, 1, 5))
})

# parameters
Theta <- mnix_sim(nSub, lambda = lambda, Omega =  Omega,
                  nu = nu * sapply(X, function(x) x[1,2]), tau = tau)

# responses
# mean vectors
Mu <- lapply(1:nSub, function(ii) {
  Theta$beta[ii,1] + Theta$beta[ii,2] * exp(-gamma * X[[ii]][,1])
})
# observation vectors
y <- lapply(1:nSub, function(ii) {
  rnorm(N[ii], mean = Mu[[ii]], sd = Theta$sigma[ii])
})
# convert data to regular format
X <- do.call(rbind, X)
y <- do.call(c, y)
id <- rep(1:nSub, times = N)

# plot data
par(mfrow = c(1,1), mar = c(4,4,.5,.5))
clrs <- rep(c("black", "blue", "red", "orange", "green4", "brown"),
            len = nSub)
plot(0, type = "n", ylab = expression(y[it]), xlab = expression(t),
     xlim = c(0, max(N)), ylim = range(sapply(Mu, range), y))
invisible(sapply(1:nSub, function(ii) {
  lines(x = 1:N[ii], y = Mu[[ii]], col = clrs[ii], lwd = 2)
  points(x = 1:N[ii], y = y[id == ii], pch = 16, cex = .8, col = clrs[ii])
}))

## ----ex2_help-----------------------------------------------------------------
# simulate hyperparameters on the transformed scale
sim_psi <- function() {
  gamma <- runif(1)
  lambda <- rnorm(2)
  Omega <- crossprod(matrix(rnorm(4), 2, 2))
  nu <- runif(1, 1, 2)
  tau <- rexp(1)/5
  list(gamma = gamma, lambda = lambda,
       logC_Omega = log_chol(Omega),
       log_nu = log(nu), log_tau = log(tau))
}
# _negative_ marginal log-posterior for ModelExt:
# R implementation
mxt_r <- function(psi, id, y, X) {
  nSub <- length(unique(id))
  lm <- sapply(1:nSub, function(ii) {
    # individual covariate and responses
    Xi <- X[id == ii,]
    yi <- y[id == ii]
    # hyperparameters on the regular scale
    gamma <- psi$gamma
    lambda <- psi$lambda
    Omega <- ilog_chol(psi$logC_Omega)
    nu <- exp(psi$log_nu) * Xi[1,2] # covariate-dependent
    tau <- exp(psi$log_tau)
    # covariate matrix
    Xi <- cbind(1, exp(-gamma * Xi[,1])) # hyperparameter-dependent
    # posterior mNIX hyperparameters
    phi_post <- mnix_post(y = yi, X = Xi,
                          lambda = lambda, Omega = Omega, nu = nu, tau = tau)
    # marginal likelihood per individual
    mnix_zeta(Omega = phi_post$Omega, nu = phi_post$nu, tau = phi_post$tau) -
      mnix_zeta(Omega = Omega, nu = nu, tau = tau)
  })
  # default prior
  lpi <- lchol_prior(ilog_chol(psi$logC_Omega))
  -(sum(lm) + lpi) # negative loglikelihood
}

## ----ex2_inst_setup, include = FALSE------------------------------------------
model_name <- "ModelExt"
p <- 2
mxt_tmb <- TMB::MakeADFun(data = c(list(model = model_name),
                                   format_data(id = id, X = X, y = y)),
                          parameters = list(gamma = 1, lambda = rep(0,p),
                                            logC_Omega = log_chol(diag(p)),
                                            log_nu = 0, log_tau = 0),
                          silent = TRUE, DLL = "losmix_TMBExports")

## ----ex2_inst, eval = FALSE---------------------------------------------------
#  # _negative marginal logposterior for ModelExt:
#  # TMB implementation
#  mxt_data <- format_data(id = id, y = y, X = X) # TMB format
#  mxt_pars <- sim_psi() # initialize with arbitrary values (placeholders)
#  mxt_tmb <- TMB::MakeADFun(data = mxt_data,
#                            parameters = mxt_pars,
#                            DLL = model_name, silent = TRUE)

## ----ex2_test-----------------------------------------------------------------
replicate(10, expr = {
  psi <- sim_psi()
  mxt_r(psi = psi, id = id, y = y, X = X) - mxt_tmb$fn(unlist(psi))
})

## ----ex2_mpost----------------------------------------------------------------
# objective function
ofun <- function(par) {
  out <- mxt_tmb$fn(par)
  attr(out, "gradient") <- mxt_tmb$gr()
  out
}

# optimization
opt <- nlm(p = mxt_tmb$par, f = ofun)
opt$code # code == 1 means that 'nlm' thinks it converged

# approximate posterior inference
npost <- 1e4 # number of posterior draws
psi_mean <- opt$estimate # (approximate) posterior mean
psi_var <- solveV(mxt_tmb$he(opt$estimate)) # (approximate) posterior variance
Psi_post <- rmvn(n = npost, mu = psi_mean, Sigma = psi_var)

## ---- eval = FALSE------------------------------------------------------------
#  mxt_tmb$simulate(psi)

## ----ex2_rxinst_setup, include = FALSE----------------------------------------
iSub <- sample(nSub, 1)
mxt1_tmb  <- TMB::MakeADFun(data = c(list(model = model_name),
                                     format_data(X = X[id == iSub,],
                                                 y = y[id == iSub])),
                            parameters = list(gamma = 1, lambda = rep(0,p),
                                              logC_Omega = log_chol(diag(p)),
                                              log_nu = 0, log_tau = 0),
                            silent = TRUE, DLL = "losmix_TMBExports")

## ----ex2_rxinst, eval = FALSE-------------------------------------------------
#  iSub <- sample(nSub, 1) # pick an observation at random
#  # TMB implementation for a single subject
#  mxt1_data <- format_data(y = y[id == iSub], X = X[id == iSub,]) # TMB format
#  mxt1_tmb <- TMB::MakeADFun(data = mxt1_data,
#                             parameters = sim_psi(), # placeholder
#                             DLL = model_name, silent = TRUE)
#  

## ----ex2_rxpost, fig.width = 7, fig.height = 3.5------------------------------
# simulate from random-effects distribution
# note: there is no need to convert to original parametrization first
Thetai_post <- apply(Psi_post, 1, function(psi) mxt1_tmb$simulate(psi))
# convert from list of lists to list with beta and sigma
Thetai_post <- unlist_bind(Thetai_post,
                           name = c("beta", "sigma"), bind = c(cbind, c))
Thetai_post$beta <- t(Thetai_post$beta) # beta needs to be transposed
thetai_true <- c(Theta$beta[iSub,], Theta$sigma[iSub]) # true parameter values

# plot
Thetai_plot <- cbind(Thetai_post$beta, Thetai_post$sigma)
thetai_names <- c("beta[i1]", "beta[i2]", "sigma[i]")
par(mfrow = c(1,3), mar = c(2,2,4,.5))
for(ii in 1:ncol(Thetai_plot)) {
  # approximate posterior
  hist(Thetai_plot[,ii], breaks = 40, xlab = "", ylab = "",
       main = parse(text = paste0("hat(p)(", thetai_names[ii],
                                  "*\" | \"*bold(Y),bold(X))")))
  # true parameter value
  abline(v = thetai_true[ii], col = "red", lwd = 2)
}

