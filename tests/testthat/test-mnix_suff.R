source("losmix-testfunctions.R")

context("mnix_suff")

test_that("Sufficient statistics are the same in R and TMB", {
  ntests <- 100
  suff_names <- c("yy", "Xy", "XX", "N")
  for(ii in 1:ntests) {
    M <- sample(1:5,1)
    N <- sample(10:20, M, replace = TRUE)
    p <- sample(1:5, 1)
    X <- sim_X(N = sum(N), p = p)
    y <- sim_y(N = sum(N))
    id <- sample(rep(sample(1000,M), times = N))
    suff_r <- tapply(1:sum(N), id,
                   function(ind) get_suff(y = y[ind], X = X[ind,]))
    suff_r <- unlist_bind(x = suff_r,
                     name = c("yy", "Xy", "XX"),
                     bind = list(c, cbind, cbind))
    suff_tmb <- mnix_suff(id = id, y = y, X = X)
    expect_equal(suff_r, suff_tmb)
  }
})
