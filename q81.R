source("common.R")

set.seed(10077)

black.scholes <- function() {
  d1 <- (log(S0/E) + (r + sigma ^ 2 / 2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  S0 * pnorm(d1) - E * exp(-r * T) * pnorm(d2)
}

grow.range <- function(r, k) {
  (r - r[1]) * k + r[1]
}

S0 <- 100
E <- 100
r <- .05
sigma <- .2
T <- .5
n <- 100

muf <- function(x,t) r * x
sigmaf <- function(x,t) sigma * x

proc.gen <- stoch.process(mu=muf, sigma=sigmaf, T)

num.simulations <- c(1000, 5000, 10000, 50000, 100000, 500000)
# num.simulations <- c(1000, 5000, 10000, 50000)
sim.price <- sapply(num.simulations, function(M) {
  S <- proc.gen(init=S0, n, M)
  
  ST <- S[, n + 1]
  ST2 <- S0 * exp(rnorm(n=M, mean=(r - sigma ^ 2 / 2) * T, sd=sigma * sqrt(T)))
  c(exp(-r * T) * c(mean(pmax(ST - E, 0)), mean(pmax(ST2 - E, 0))), black.scholes())
})

pdf("figure11.pdf")
x <- 1:length(num.simulations)
plot(x, sim.price[1,], type='b', ylim=grow.range(range(sim.price), 1.2), 
     col="darkblue", lwd=2,
     pch=17, xaxt='n', xlab="Number of simulations", ylab="Option price")
lines(x, sim.price[2,], type='b', col="darkgreen", lwd=2, pch=18)
lines(x, sim.price[3,], type='b', col="red", lwd=2, pch=15)
axis(1, x, num.simulations)
legend('topright', c("Discretisation", "Formula 6.16", "Black-Scholes"),
       col=c("darkblue", "darkgreen", "red"), lwd=2)
dev.off()
d <- abs(sim.price - black.scholes())
