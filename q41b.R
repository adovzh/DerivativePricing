source("common.R")

mu <- .15
sigma <- .2
T <- 2
n <- 100
M=1000

Y <- stoch.process.const(mu=mu-sigma*sigma/2, sigma=sigma, T=T)(init=0, n=n, M=M)
X <- exp(Y)
x <- seq(from=0, to=T, by=T/n)
cols <- rainbow(10, alpha=.8)
pdf("figure04.pdf", width=14, heigh=7)
p <- par(no.readonly=TRUE)
par(mfrow=c(1, 2))

plot(x, X[1,], type="l", 
     main="Sample paths for Geometric Brownian Motion process",
     xlab="time t", ylab="x(t)", ylim=range(X))
for (i in 2:M) lines(x, X[i,], col=cols[i %% length(cols)])

plot(x, Y[1,], type="l", 
     main="Sample paths for a log return process",
     xlab="time t", ylab="process", ylim=range(Y))
for (i in 2:M) lines(x, Y[i,], col=cols[i %% length(cols)])
par(p)
dev.off()

# ii) 
pdf("figure05.pdf", width=14, height=7)
p <- par(no.readonly=TRUE)
par(mfrow=c(1, 2))
hist(X[, n+1], breaks=50, density=20, prob=TRUE,
     main="Sample distribution at t=2", xlab="Log GBM process")
curve(dlnorm(x, meanlog=(mu-sigma*sigma/2)*T, sdlog=sigma*sqrt(T)), col="darkblue", lwd=2, add=TRUE)
hist(Y[, n+1], breaks=50, density=20, prob=TRUE, 
     main="Sample distribution at t=1", xlab="GBM process")
curve(dnorm(x, mean=(mu-sigma*sigma/2)*T, sd=sigma*sqrt(T)), col="darkblue", lwd=2, add=TRUE)
par(p)
dev.off()

# iii)
lgbm <- stoch.process.const(mu=mu-sigma*sigma/2, sigma=sigma, T=T)
num.intervals <- c(100, 1000, 5000)
num.simulations <- c(1000, 10000, 50000)
pdf("figure06.pdf")
p <- par(no.readonly=TRUE)
par(mfrow=c(length(num.intervals), length(num.simulations)))

intervals <- rep(num.intervals, each=length(num.simulations))
simulations <- rep(num.simulations, times=length(num.intervals))
r <- mapply(function(i,j) {
  Y <- lgbm(init=0, n=i, M=j)
  X <- exp(Y)
  hist(X[, i + 1], breaks=50, density=20, prob=TRUE,
       main=paste("n =",i,"; M =", j), xlab="GBM process")
  curve(dlnorm(x, meanlog=(mu-sigma*sigma/2)*T, sdlog=sigma*sqrt(T)), 
        col="darkblue", lwd=2, add=TRUE)    
  list("name"=sprintf("n=%d, M=%d", i, j), 
       "mean of X(2)"=mean(X[, i+1]),
       "mean of Y(2)"=mean(Y[,i+1]),
       "variance of X(2)"=var(X[, i+1]),
       "variance of Y(2)"=var(Y[, i+1]))
}, intervals, simulations)
r["name",] -> colnames(r)
r <- r[rownames(r) != "name", ]
r <- cbind(c(exp(mu*T), 
                         (mu-sigma^2/2)*T, 
                         exp(2*mu*T)*(exp(sigma^2*T)-1),
                         sigma^2*T), r)
par(p)
dev.off()
