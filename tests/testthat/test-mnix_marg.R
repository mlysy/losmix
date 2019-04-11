source("losmix-testfunctions.R")

context("mnix_marg")

test_that("R(loglik + lprior) = TMB(lmarg + lcond)", {
  ntests <- 20
  loglik <- function(theta, y, X) {
    sum(dnorm(y, mean = X %*% theta$beta, sd = sqrt(theta$sig2), log = TRUE))
  }
  logpi <- function(theta, Phi) {
    dmnix(x = theta$beta, v = theta$sig2, Phi = Phi, log = TRUE)
  }
  ans <- rep(NA, ntests)
  for(ii in 1:ntests) {
    M <- sample(1:5,1)
    N <- sample(10:20, M, replace = TRUE)
    p <- sample(1:5, 1)
    X <- lapply(N, sim_X, p = p)
    y <- lapply(N, sim_y)
    theta <- replicate(M, sim_theta(p = p), simplify = FALSE)
    Phi <- sim_Phi(p = p)
    # R calculation
    suff_r <- mapply(get_suff, y = y, X = X, SIMPLIFY = FALSE)
    ll <- sum(mapply(loglik, y = y, X = X, theta = theta))
    lpi <- sum(sapply(theta, logpi, Phi = Phi))
    lpi <-  lpi + lchol_prior(Phi$Omega) # default prior on Phi
    # TMB calculation
    id <- sim_id(N)
    y <- do.call(c, y)[id$order]
    X <- do.call(rbind, X)[id$order,]
    id <- id$id
    post_tmb <- do.call(mnix_post, c(list(id = id, y = y, X = X), Phi))
    lc <- sapply(1:M, function(ii) {
      logpi(theta = theta[[ii]], Phi = list(lambda = post_tmb$lambda[ii,],
                                            Omega = post_tmb$Omega[,,ii],
                                            tau = post_tmb$tau[ii],
                                            nu = post_tmb$nu[ii]))
    })
    lc <- sum(lc)
    omarg <- mnix_marg(id = id, y = y, X = X)
    # TMB returns NLL up to factor of N/2 * log(2*pi)
    lm <- -omarg$fn(vec_phi(Phi)) - sum(N/2 * log(2*pi))
    expect_equal(ll + lpi, lc + lm)
  }
})


#--- scratch ---------------------------------------------------------------

## export_compile <- function() {
##   RcppTMBTest::export_models()
##   pkgbuild::compile_dll(force = TRUE)
##   devtools::install(quick = TRUE)
## }
## export_compile()

## q()

## require(losmix)
## require(testthat)
## source("losmix-testfunctions.R")

## lchol_prior2 <- function(logC) {
##   p <- (-1 + sqrt(1 + 8*length(logC)))/2
##   jj <- 0
##   lpi <- rep(NA, p-1)
##   for(ii in 1:(p-1)) {
##     lpi[ii] <- logC[jj+1]
##     jj <- jj+(ii+1)
##   }
##   lpi
## }

## M <- sample(1:5,1)
## N <- sample(10:20, M, replace = TRUE)
## p <- sample(1:5, 1)
## X <- lapply(N, sim_X, p = p)
## y <- lapply(N, sim_y)
## id <- sim_id(N)
## y <- do.call(c, y)[id$order]
## X <- do.call(rbind, X)[id$order,]
## id <- id$id
## omarg <- mnix_marg(id = id, y = y, X = X)
## Phi <- sim_Phi(p = p)
## omarg$report(vec_phi(Phi))$zeta1 - omarg$report(vec_phi(Phi))$zeta2
