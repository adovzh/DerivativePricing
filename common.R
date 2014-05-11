stoch.process <- function(mu, sigma, T) {
  function(init, n, M) {
    dt <- T / n
    X <- rep(init, M)
    dim(X) <- c(M, 1)
    A <- X
    t <- 0
    
    for (i in seq(1, T * n)) {
      dx <- mapply(mu, X, rep(t, M)) * dt + 
        mapply(sigma, X, rep(t, M)) * rnorm(M, sd=sqrt(dt))
      X <- X + dx
      t <- t + dt
      A <- cbind(A, X)
    }
    
    A
  }
}

# mu & sigma are constants here
stoch.process.const <- function(mu, sigma, T) {
  function(init, n, M) {
    dt <- T / n
    A <- rnorm(n * M, mean=mu * dt, sd=sigma * sqrt(dt))
    dim(A) <- c(M, n)
    A <- cbind(rep(init, M), A)
    A <- t(apply(A, 1, cumsum))
  }
}