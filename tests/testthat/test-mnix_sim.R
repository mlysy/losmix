## require(losmix)
## require(testthat)
source("losmix-testfunctions.R")


context("mnix_sim")

ntests <- 20

test_that("R == TMB: arbitrary mNIX, single params", {
  seeds <- sample(1000,ntests)
  for(ii in 1:ntests) {
    p <- sample(1:5, 1)
    n <- sample(1:5, 1)
    Phi <- sim_Phi(p)
    # R calculation
    set.seed(seeds[ii])
    theta_r <- replicate(n, do.call(rmnix, Phi), simplify = FALSE)
    theta_r <- unlist_bind(theta_r, c("beta", "sigma"), c(rbind, c))
    # TMB calculation
    set.seed(seeds[ii])
    theta_tmb <- do.call(mnix_sim, c(list(n = n), Phi))
    expect_equal(theta_tmb, theta_r)
  }
})


test_that("R == TMB: arbitrary mNIX, multiple params", {
  for(ii in 1:ntests) {
    p <- sample(1:5, 1)
    n <- sample(1:5, 1)
    Phi <- replicate(n, sim_Phi(p), simplify = FALSE)
    # R calculation
    set.seed(1)
    theta_r <- lapply(Phi, function(phi) do.call(rmnix, phi))
    theta_r <- unlist_bind(theta_r, c("beta", "sigma"), c(rbind, c))
    # TMB calculation
    set.seed(1)
    Phi_tmb <- unlist_bind(Phi, c("lambda", "Omega", "nu", "tau"),
                           c(rbind, cbind, c, c))
    Phi_tmb$Omega <- array(Phi_tmb$Omega, c(p,p,n))
    theta_tmb <- do.call(mnix_sim, c(list(n = n), Phi_tmb))
    expect_equal(theta_tmb, theta_r)
  }
})

test_that("R == TMB: posterior mNIX, single params", {
  seeds <- sample(1000,ntests)
  for(ii in 1:ntests) {
    N <- sample(10:20, 1)
    p <- sample(1:5, 1)
    X <- sim_X(N = N, p = p)
    y <- sim_y(N = N)
    n <- sample(1:5, 1)
    Phi <- sim_Phi(p)
    # R calculation
    set.seed(seeds[ii])
    Phi_post <- get_post(suff = get_suff(y = y, X = X), Phi = Phi)
    theta_r <- replicate(n, do.call(rmnix, Phi_post), simplify = FALSE)
    theta_r <- unlist_bind(theta_r, c("beta", "sigma"), c(rbind, c))
    # TMB calculation
    set.seed(seeds[ii])
    theta_tmb <- do.call(mnix_sim, c(list(n = n, y = y, X = X), Phi))
    expect_equal(theta_tmb, theta_r)
  }
})

test_that("R == TMB: posterior mNIX, multiple params", {
  seeds <- sample(1000,ntests)
  for(ii in 1:ntests) {
    N <- sample(10:20, 1)
    p <- sample(1:5, 1)
    X <- sim_X(N = N, p = p)
    y <- sim_y(N = N)
    n <- sample(1:5, 1)
    Phi <- replicate(n, sim_Phi(p), simplify = FALSE)
    # R calculation
    set.seed(seeds[ii])
    suff <- get_suff(y = y, X = X)
    Phi_post <- lapply(Phi, get_post, suff = suff)
    theta_r <- lapply(Phi_post, function(phi) do.call(rmnix, phi))
    theta_r <- unlist_bind(theta_r, c("beta", "sigma"), c(rbind, c))
    # TMB calculation
    set.seed(seeds[ii])
    Phi_tmb <- unlist_bind(Phi, c("lambda", "Omega", "nu", "tau"),
                           c(rbind, cbind, c, c))
    Phi_tmb$Omega <- array(Phi_tmb$Omega, c(p,p,n))
    theta_tmb <- do.call(mnix_sim, c(list(n = n, y = y, X = X), Phi_tmb))
    expect_equal(theta_tmb, theta_r)
  }
})


test_that("mnix_sim and mnix_marg$simulate give identical resuls", {
  seeds <- sample(1000,ntests)
  for(ii in 1:ntests) {
    M <- sample(1:5, 1)
    N <- sample(10:20, M, replace = TRUE)
    p <- sample(1:5, 1)
    X <- do.call(rbind, lapply(N, sim_X, p = p))
    y <- do.call(c, lapply(N, sim_y))
    id <- rep(1:M, times = N)
    n <- sample(1:5, 1)
    Phi <- replicate(n, sim_Phi(p), simplify = FALSE)
    # mnix_sim calculation
    set.seed(seeds[ii])
    Theta_sim <- lapply(Phi, function(phi) {
      do.call(mnix_sim, c(list(n = M, y = y, X = X, id = id), phi))
    })
    # mnix_marg calculation
    set.seed(seeds[ii])
    nll <- mnix_marg(id = id, y = y, X = X)
    Theta_marg <- apply(sapply(Phi, vec_phi), 2, function(psi) {
      ans <- nll$simulate(psi)
      list(beta = t(ans$beta), sigma = ans$sigma)
    })
    expect_equal(Theta_sim, Theta_marg)
  }
})
