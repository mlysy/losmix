# generate a list of lists
Phi <- replicate(5, {
  list(lambda = rnorm(3), Omega = diag(rnorm(3)),
       nu = rexp(1), tau = rexp(1))
}, simplify = FALSE)

unlist_bind(Phi, bind = c(cbind, rbind, c, c))
