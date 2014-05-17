source("common.R")

set.seed(1004)

k <- 0.5
equilibrium <- 0.065
mu <- function(x,t) k * (equilibrium - x)
sigma <- 0.02
sigmaf <- function(x,t) sigma
T <- 1

x0 <- 0.06
n <- 100
M <- 1000

ou <- stoch.process(mu=mu, sigma=sigmaf, T=T)(init=x0, n=n, M=M)

x <- seq(from=0, to=T, by=T/n)
cols <- rainbow(10, alpha=.8)
pdf("figure07.pdf")
plot(x, ou[1,], type="l", 
     main="Sample paths for a O-U process",
     xlab="time t", ylab="O-U process", ylim=range(ou))
for (i in 2:M) lines(x, ou[i,], col=cols[i %% length(cols)])
dev.off()

ou.mean <- ou.mean.gen(k, equilibrium, x0)
ou.sd <- ou.sd.gen(k, sigma)

pdf("figure08.pdf", width=14, height=7)
p <- par(no.readonly=TRUE)
par(mfrow=c(1, 2))
hist(ou[, n/2+1], breaks=50, density=20, prob=TRUE,
     main="Sample distribution at t=6 months", xlab="O-U process")
curve(dnorm(x, 
            mean=ou.mean(T/2),
            sd=ou.sd(T/2)
            ), col="darkblue", lwd=2, add=TRUE)
hist(ou[, n+1], breaks=50, density=20, prob=TRUE, 
     main="Sample distribution at t=12 months", xlab="O-U process")
curve(dnorm(x, 
            mean=ou.mean(T),
            sd=ou.sd(T)
            ), 
      col="darkblue", lwd=2, add=TRUE)
par(p)
dev.off()

r <- data.frame(theoretical=c(ou.mean(0.5), ou.mean(1),
                              ou.sd(0.5)^2, ou.sd(1)^2),
                actual=c(mean(ou[, n/2+1]), mean(ou[, n+1]),
                         var(ou[, n/2+1]), var(ou[, n+1])),
                names=c("mean of OU(1/2)", 
                        "mean of OU(1)", 
                        "variance of OU(1/2)",
                        "variance of OU(1)"),
                row.names="names")
write.csv(r, "r3.csv")
