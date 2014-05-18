set.seed(1005)

ynt.gen <- function(k, n) {
  function(t) {
    dt <- t / n
    dfc <- sapply(seq(0, n - 1), function(i) exp(-k * (t - i * dt)))
    sum(dfc * rnorm(n, mean=0, sd=sqrt(dt)))
  }
}

theor.var <- function(k, t) {
  (1 - exp(-2 * k * t)) / (2 * k)
}

k <- .5
t <- 1
n <- 100
M <- 1000

Yn <- ynt.gen(k, n)
Y <- replicate(M, Yn(t))

tmean <- 0
tvar <- theor.var(k, t)

pdf("figure09.pdf")
hist(Y, breaks=50, density=20, probability=TRUE,
     main="Stochastic integral", xlab="time t")
curve(dnorm(x, mean=tmean, sd=sqrt(tvar)),
      col="darkblue", lwd=2, add=TRUE)

num.intervals <- c(100, 5000)
num.simulations <- c(1000, 10000)
dev.off()

pdf("figure10.pdf")
p <- par(no.readonly=TRUE)
par(mfrow=c(length(num.intervals), length(num.simulations)))

intervals <- rep(num.intervals, each=length(num.simulations))
simulations <- rep(num.simulations, times=length(num.intervals))

report <- mapply(function(n,M) {
  Yn <- ynt.gen(k, n)
  Y <- replicate(n=M, expr=Yn(t))
  
  hist(Y, breaks=50, density=20, probability=TRUE,
       main=sprintf("Yn(1) (n=%d, M=%d)", n, M), 
       xlab="time t")
  curve(dnorm(x, mean=tmean, sd=sqrt(tvar)),
        col="darkblue", lwd=2, add=TRUE)
  list("names"=sprintf("n=%d, M=%d", n, M),
       "Mean of Y(1)"=mean(Y),
       "Variance of Y(1)"=var(Y))
}, intervals, simulations)

par(p)
dev.off()

report["names",] -> colnames(report)
report <- report[rownames(report) != "names", ]
report <- cbind(Theoretical=c(tmean, tvar), report)
                

# report <- data.frame(Theoretical=c(0, theor.var(k, t)),
#                      Actual=c(mean(Y), var(Y)),
#                      names=c("Mean of Y(1)", "Variance of Y(1)"),
#                      row.names="names")
write.csv(report, file="r5.csv")
