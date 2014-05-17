# mu and sigma are functions (x,t) here
stoch.process <- function(mu, sigma, T) {
  function(init, n, M) {
    dt <- T / n
    X <- rep(init, M)
    dim(X) <- c(M, 1)
    A <- X
    t <- 0
    
    for (i in seq(1, n)) {
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

# theoretical mean function of O-U process
ou.mean.gen <- function(k, mu, x0) {
  function(t) {
    x0*exp(-k*t) + mu*(1-exp(-k*t))
  }
}

# theoretical standard deviation of O-U process
ou.sd.gen <- function(k, sigma) {
  function(t) {
    sqrt(sigma^2/(2*k)*(1-exp(-2*k*t)))
  }
}