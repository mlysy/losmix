## require(losmix)
## require(testthat)
source("losmix-testfunctions.R")

context("mnix_post")

test_that("R == TMB: single data, multiple params", {
  ntests <- 100
  for(ii in 1:ntests) {
    M <- sample(1:5,1)
    N <- sample(10:20, M, replace = TRUE)
    p <- sample(1:5, 1)
    X <- sim_X(N = sum(N), p = p)
    y <- sim_y(N = sum(N))
    id <- rep(1, sum(N))
    Phi <- replicate(M, sim_Phi(p = p), simplify = FALSE)
    # R calculation
    suff_r <- get_suff(y = y, X = X)
    post_r <- lapply(Phi, get_post, suff = suff_r)
    post_r <- unlist_bind(x = post_r,
                          name = c("lambda", "Omega", "tau", "nu"),
                          bind = list(rbind, cbind, c, c))
    post_r$Omega <- array(post_r$Omega, dim = c(p, p, M))
    # TMB calculation
    Phi_tmb <- unlist_bind(x = Phi,
                           name = c("lambda", "Omega", "tau", "nu"),
                           bind = list(rbind, cbind, c, c))
    Phi_tmb$Omega <- array(Phi_tmb$Omega, dim = c(p, p, M))
    post_tmb <- do.call(mnix_post, c(list(id = id, y = y, X = X), Phi_tmb))
    expect_equal(post_r, post_tmb)
  }
})

test_that("R == TMB: multiple data, single params", {
  ntests <- 100
  for(ii in 1:ntests) {
    M <- sample(1:5,1)
    N <- sample(10:20, M, replace = TRUE)
    p <- sample(1:5, 1)
    X <- sim_X(N = sum(N), p = p)
    y <- sim_y(N = sum(N))
    id <- sample(rep(sample(1000,M), times = N))
    Phi <- sim_Phi(p)
    # R calculation
    suff_r <- tapply(1:sum(N), id,
                     function(ind) get_suff(y = y[ind], X = X[ind,]))
    post_r <- lapply(suff_r, get_post, Phi = Phi)
    post_r <- unlist_bind(x = post_r,
                          name = c("lambda", "Omega", "tau", "nu"),
                          bind = list(rbind, cbind, c, c))
    post_r$Omega <- array(post_r$Omega, dim = c(p, p, M))
    # TMB calculation
    post_tmb <- do.call(mnix_post, c(list(id = id, y = y, X = X), Phi))
    expect_equal(post_r, post_tmb)
  }
})

test_that("R == TMB: multiple data and params", {
  ntests <- 100
  for(ii in 1:ntests) {
    M <- sample(1:5,1)
    N <- sample(10:20, M, replace = TRUE)
    p <- sample(1:5, 1)
    X <- sim_X(N = sum(N), p = p)
    y <- sim_y(N = sum(N))
    id <- sample(rep(sample(1000,M), times = N))
    Phi <- replicate(M, sim_Phi(p = p), simplify = FALSE)
    # R calculation
    suff_r <- tapply(1:sum(N), id,
                     function(ind) get_suff(y = y[ind], X = X[ind,]))
    post_r <- mapply(get_post, suff = suff_r, Phi = Phi, SIMPLIFY = FALSE)
    post_r <- unlist_bind(x = post_r,
                          name = c("lambda", "Omega", "tau", "nu"),
                          bind = list(rbind, cbind, c, c))
    post_r$Omega <- array(post_r$Omega, dim = c(p, p, M))
    # TMB calculation
    Phi_tmb <- unlist_bind(x = Phi,
                           name = c("lambda", "Omega", "tau", "nu"),
                           bind = list(rbind, cbind, c, c))
    Phi_tmb$Omega <- array(Phi_tmb$Omega, dim = c(p, p, M))
    post_tmb <- do.call(mnix_post, c(list(id = id, y = y, X = X), Phi_tmb))
    expect_equal(post_r, post_tmb)
  }
})
