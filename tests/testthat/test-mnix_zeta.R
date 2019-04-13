require(losmix)
require(testthat)
source("losmix-testfunctions.R")

context("mnix_zeta")

test_that("mnix_zeta and lchol_prior are same in TMB and R", {
  ntests <- 20
  for(ii in 1:ntests) {
    M <- sample(1:5,1)
    N <- sample(10:20, M, replace = TRUE)
    p <- sample(1:5, 1)
    X <- do.call(rbind, lapply(N, sim_X, p = p))
    y <- do.call(c, lapply(N, sim_y))
    id <- rep(1:M, times = N)
    Phi <- sim_Phi(p = p)
    # R calculation
    Phi_post <- do.call(mnix_post, c(list(y = y, X = X, id = id), Phi))
    lm_r <- do.call(mnix_zeta, Phi[c("Omega", "nu", "tau")])
    lm_r <- do.call(mnix_zeta, Phi_post[c("Omega", "nu", "tau")])
    lm_r <- c(lm_r, do.call(mnix_zeta, Phi[c("Omega", "nu", "tau")]))
    lm_r <- c(lm_r, lchol_prior(Phi$Omega))
    lm_r <- do.call(mnix_zeta, Phi_post[c("Omega", "nu", "tau")])
    lm_r <- sum(lm_r - do.call(mnix_zeta, Phi[c("Omega", "nu", "tau")]))
    lm_r <- lm_r + lchol_prior(Phi$Omega)
    # TMB calculation
    nll_tmb <- mnix_marg(id = id, y = y, X = X)
    lm_tmb <- -nll_tmb$fn(vec_phi(Phi))
    expect_equal(lm_r, lm_tmb)
  }
})

