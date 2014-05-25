set.seed(1003)

mu <- .15
sigma <- .2
T <- 2

num.simulations <- c(1000, 10000, 50000)

pdf("figure12.pdf", width=10.5)
p <- par(no.readonly=TRUE)
par(mfcol=c(2, length(num.simulations)))

r2 <- sapply(num.simulations, function(M) {
  y <- rnorm(n=M, mean=(mu - sigma ^ 2 / 2) * T, sd=sigma * sqrt(T))
  x <- exp(y)
  hist(x, breaks=50, density=20, prob=TRUE,
       main=paste(sprintf("X(2) M = %d", M)), xlab="X")
  curve(dlnorm(x, meanlog=(mu-sigma*sigma/2)*T, sdlog=sigma*sqrt(T)), 
        col="darkblue", lwd=2, add=TRUE)    
  
  hist(y, breaks=50, density=20, prob=TRUE,
       main=paste(sprintf("Y(2) M = %d", M)), xlab="Y")
  curve(dnorm(x, mean=(mu-sigma*sigma/2)*T, sd=sigma*sqrt(T)), 
        col="darkblue", lwd=2, add=TRUE)    
  
  list("name"=sprintf("M=%d", M), 
       "mean of X(2)"=mean(x),
       "mean of Y(2)"=mean(y),
       "variance of X(2)"=var(x),
       "variance of Y(2)"=var(y))
})

par(p)
dev.off()

r2["name", ] -> colnames(r2)
r2 <- r2[rownames(r2) != "name", ]
r2 <- cbind(Theoretical=c(exp(mu*T), 
                         (mu-sigma^2/2)*T, 
                         exp(2*mu*T)*(exp(sigma^2*T)-1),
                         sigma^2*T), r2)
