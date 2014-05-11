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

set.seed(1002)
wiener <- stoch.process(mu = function(x,t) 0, 
                        sigma = function(x,t) 1,
                        T = 1)
P1 <- wiener(init = 1, n = 100, M = 1000)

x <- seq(0, 1, 1 / 100)
cols <- rainbow(10, alpha=.8)
plot(x, P1[1,], type="l", 
     main="Sample paths for a pure Wiener process",
     xlab="time t", ylab="Wiener process", ylim=range(P1))
for (i in 2:100) lines(x, P1[i,], col=cols[i %% length(cols)])

hist(P1[,51], breaks=50, density=20, prob=T,
     main="Sample distribution at t=0.5", xlab="Wiener process")
curve(dnorm(x, mean=1, sd=sqrt(.5)), col="darkblue", lwd=2, add=T)
hist(P1[,101], breaks=50, density=20, prob=T, 
     main="Sample distribution at t=1", xlab="Wiener process")
curve(dnorm(x, mean=1, sd=1), col="darkblue", lwd=2, add=T)

# ii)
p <- par(no.readonly=T)
par(mfrow=c(2,2))
for (i in c(100, 1000)) {
  for (j in c(1000, 10000)) {
    P <- wiener(init = 1, n = i, M = j)
    hist(P[, i / 2 + 1], breaks=50, density=20, prob=T,
         main=paste("n =",i,"; M =", j), xlab="Wiener process")
    curve(dnorm(x, mean=1, sd=sqrt(.5)), col="darkblue", lwd=2, add=T)
  }
}
par(p)
