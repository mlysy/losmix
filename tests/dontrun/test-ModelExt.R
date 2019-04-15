require(losmix)
require(TMB)

model_name <- "ModelExt"
include_path <- system.file("include", package = "losmix")
TMB::compile(paste0(model_name, ".cpp"),
             PKG_CXXFLAGS = paste0('-I"', include_path, '"'))
dyn.load(TMB::dynlib(model_name))

# simulate data

nSub <- 10 # number of subjects
N <- sample(20:50, nSub, replace = TRUE) # observations per subject

# hyperparameters
gamma <- runif(1, .1, .5)
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
par(mfrow = c(1,1))
plot(0, type = "n", ylab = "y", xlab = "t",
     xlim = c(0, max(N)), ylim = range(sapply(Mu, range), y))
invisible(sapply(1:nSub, function(ii) {
  lines(x = 1:N[ii], y = Mu[[ii]])
  points(x = 1:N[ii], y = y[id == ii], pch = 16, cex = .5)
}))


# check that calculations are correct

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

# _negative_ marginal loglikelihood of ModelExt in R
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


nll_tmb <- TMB::MakeADFun(data = format_data(id = id, y = y, X = X),
                          parameters = sim_psi(),
                          DLL = model_name, silent = TRUE)

replicate(10, expr = {
  psi <- sim_psi()
  mxt_marg(psi = psi, id = id, y = y, X = X) - nll_tmb$fn(unlist(psi))
})

nll_sim <- TMB::MakeADFun(data = c(list(model_name = "ModelExt_marg"),
                                   format_data(y = y[id == 1],
                                               X = X[id == 1,])),
                          parameters = sim_psi(),
                          DLL = model_name, silent = TRUE)

system.time({
  replicate(1e4, nll_sim$simulate(unlist(psi)))
})

#--- scratch -------------------------------------------------------------------

nuX_r <- exp(psi$log_nu) * tapply(X[,2], id, function(x) x[1])
nuX_tmb <- nll_tmb$report(unlist(psi))$nuX
range(nuX_r - nuX_tmb)
XtrGamma_r <- rbind(1, exp(-psi$gamma * X[,1]))
XtrGamma_tmb <- nll_tmb$report(unlist(psi))$XtrGamma
range(XtrGamma_r - XtrGamma_tmb)


psi <- sim_psi()
nll_tmb$report(unlist(psi))$lm
mxt_marg(psi = psi, id = id, y = y, X = X)
