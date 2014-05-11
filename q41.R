source("common.R")

set.seed(1002)
wiener <- stoch.process(mu = function(x,t) 0, 
                        sigma = function(x,t) 1,
                        T = 1)
P1 <- wiener(init = 1, n = 100, M = 1000)

x <- seq(0, 1, 1 / 100)
cols <- rainbow(10, alpha=.8)
pdf("figure01.pdf")
plot(x, P1[1,], type="l", 
     main="Sample paths for a pure Wiener process",
     xlab="time t", ylab="Wiener process", ylim=range(P1))
for (i in 2:1000) lines(x, P1[i,], col=cols[i %% length(cols)])
dev.off()

pdf("figure02.pdf", width=14, height=7)
p <- par(no.readonly=TRUE)
par(mfrow=c(1, 2))
hist(P1[,51], breaks=50, density=20, prob=TRUE,
     main="Sample distribution at t=0.5", xlab="Wiener process")
curve(dnorm(x, mean=1, sd=sqrt(.5)), col="darkblue", lwd=2, add=TRUE)
hist(P1[,101], breaks=50, density=20, prob=TRUE, 
     main="Sample distribution at t=1", xlab="Wiener process")
curve(dnorm(x, mean=1, sd=1), col="darkblue", lwd=2, add=TRUE)
par(p)
dev.off()

# ii)
wiener <- stoch.process.const(mu=0, sigma=1, T=1)
num.intervals <- c(100, 1000, 5000)
num.simulations <- c(1000, 10000, 50000)
pdf("figure03.pdf")
p <- par(no.readonly=TRUE)
par(mfrow=c(length(num.intervals), length(num.simulations)))

intervals <- rep(num.intervals, each=length(num.simulations))
simulations <- rep(num.simulations, times=length(num.intervals))
r <- mapply(function(i,j) {
  P <- wiener(init=1, n=i, M=j)
  hist(P[, i + 1], breaks=50, density=20, prob=TRUE,
       main=paste("n =",i,"; M =", j), xlab="Wiener process")
  curve(dnorm(x, mean=1, sd=1), col="darkblue", lwd=2, add=TRUE)    
  list("name"=sprintf("n=%d, M=%d", i, j), 
       "mean at t=0.5"=mean(P[, i/2+1]),
       "mean at t=1"=mean(P[,i+1]),
       "variance at t=0.5"=var(P[, i/2+1]),
       "variance at t=1"=var(P[, i+1]))
}, intervals, simulations)
r["name",] -> colnames(r)
r <- r[rownames(r) != "name", ]
par(p)
dev.off()