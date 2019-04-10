source("losmix-testfunctions.R")

context("get_tmbdata")

test_that("Extraction is as expected", {
  ntests <- 100
  for(ii in 1:ntests) {
    M <- sample(1:5,1)
    N <- sample(10:20, M, replace = TRUE)
    p <- sample(1:5, 1)
    X <- sim_X(N = sum(N), p = p)
    y <- sim_y(N = sum(N))
    lab <- sample(1000,M)
    id <- sample(rep(lab, times = N))
    tmblist <- losmix:::get_tmbdata(id = id, y = y, X = X)
    jj <- sample(length(tmblist$nObs), 1)
    ljj <- names(tmblist$nObs) == lab[jj]
    expect_equal(y[id == lab[jj]],
                 tmblist$y[tmblist$iStart[ljj]+1:tmblist$nObs[ljj]])
    expect_equal(X[id == lab[jj],,drop=FALSE],
                 tmblist$X[tmblist$iStart[ljj]+1:tmblist$nObs[ljj],,drop=FALSE])
  }
})
