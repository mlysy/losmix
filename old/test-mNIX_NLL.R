require(losmix)
require(testthat)
## require(numDeriv)
source("losmix-testfunctions.R")

context("Single mNIX")

test_that("Sufficient statistics are the same in R and TMB", {
  ntests <- 10
  suff_names <- c("yy", "Xy", "XX", "N")
  for(ii in 1:ntests) {
    N <- sample(10:20, 1)
    p <- sample(3:5, 1)
    X <- sim_X(N, p)
    y <- sim_y(N)
    odata <- list(Xtr = t(X), y = y)
    opars <- Phi2par(Phi = sim_Phi(p))
    obj <- TMB::MakeADFun(data = c(list(model_name = "mNIX_NLL"), odata),
                          parameters = opars,
                          DLL = "losmix_TMBExports", silent = TRUE)
    suff_r <- get_suff(y, X)
    suff_tmb <- obj$simulate(par = vec_phi(sim_Phi(p)))[suff_names]
    expect_equal(suff_r, suff_tmb)
  }
})

test_that("Conjugate posterior hyperparameters are the same in R and TMB", {
  ntests <- 10
  Phi_names <- c("lambda", "Omega", "nu", "tau")
  names(Phi_names) <- paste0(Phi_names, "_hat")
  for(ii in 1:ntests) {
    N <- sample(10:20, 1)
    p <- sample(3:5, 1)
    X <- sim_X(N, p)
    y <- sim_y(N)
    Phi <- sim_Phi(p)
    odata <- list(Xtr = t(X), y = y)
    opars <- Phi2par(sim_Phi(p))
    obj <- TMB::MakeADFun(data = c(list(model_name = "mNIX_NLL"), odata),
                          parameters = opars,
                          DLL = "losmix_TMBExports", silent = TRUE)
    Phi_hat_r <- get_post(suff = get_suff(y, X), Phi = Phi)
    Phi_hat_r$lambda <- as.matrix(Phi_hat_r$lambda)
    Phi_hat_tmb <- obj$simulate(par = vec_phi(Phi))[names(Phi_names)]
    names(Phi_hat_tmb) <- as.character(Phi_names)
    expect_equal(Phi_hat_r, Phi_hat_tmb)
  }
})

test_that("R(loglik + lprior) = TMB(lmarg + lcond)", {
  ntests <- 10
  Phi_names <- c("lambda", "Omega", "nu", "tau")
  names(Phi_names) <- paste0(Phi_names, "_hat")
  loglik <- function(theta, y, X) {
    sum(dnorm(y, mean = X %*% theta$beta, sd = sqrt(theta$sig2), log = TRUE))
  }
  logpi <- function(theta, Phi) {
    dmnix(x = theta$beta, v = theta$sig2, Phi = Phi, log = TRUE)
  }
  ans <- rep(NA, ntests)
  for(ii in 1:ntests) {
    N <- sample(10:20, 1)
    p <- sample(3:5, 1)
    X <- sim_X(N, p)
    y <- sim_y(N)
    theta <- sim_theta(p)
    Phi <- sim_Phi(p)
    odata <- list(Xtr = t(X), y = y)
    opars <- Phi2par(sim_Phi(p))
    obj <- TMB::MakeADFun(data = c(list(model_name = "mNIX_NLL"), odata),
                          parameters = opars,
                          DLL = "losmix_TMBExports", silent = TRUE)
    Phi_hat <- obj$simulate(par = vec_phi(Phi)) # conjugate posterior
    Phi_hat <- setNames(Phi_hat[names(Phi_names)], as.character(Phi_names))
    Phi_hat$lambda <- c(Phi_hat$lambda)
    ll <- loglik(theta = theta, y = y, X = X)
    lpi <- logpi(theta = theta, Phi = Phi)
    lc <- logpi(theta, Phi = Phi_hat)
    # TMB returns NLL up to factor of N/2 * log(2*pi)
    lm <- -obj$fn(vec_phi(Phi)) - N/2 * log(2*pi)
    expect_equal(ll + lpi, lc + lm)
  }
})

