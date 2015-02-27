
#1
GMM <- function() {
  X <- runif(1)
  if (X <= .3) {
    return(rnorm(1, mean = 1, sd = 1))
  } else {
    return(rnorm(1, mean = 3, sd = sqrt(1.5)))
  }
  return(X)
}

samples <- sapply(1:100, function(x) GMM())
plot(density(samples), xlab='Kernel Density plot of X',main='')

#2

likelihood <- function(p,x = samples) {
  L <- prod((p*dnorm(x, mean = 1, sd = 1) + (1-p)*dnorm(x, mean = 3, sd = sqrt(1.5))))
  return(L)
}

p <- seq(0,1,.001)
samples2 <- sapply(p, function(p) likelihood(p=p))
qplot(p, samples2, main = "Full Path Over Pi") + xlab("Values of pi") + ylab("Likelihood")

#3

g.pi <- 1
x <- samples

integrand <- Vectorize(function(p) {1 * prod(p*dnorm(x, mean = 1, sd = 1) + (1-p)*dnorm(x, mean = 3, sd = sqrt(1.5)))})
integrate(integrand, lower=0, upper=1)

#4

x <- samples

g <- expression(6*x^5)
product <- expression(prod(p*dnorm(x, mean = 1, sd = 1) + (1-p)*dnorm(x, mean = 3, sd = sqrt(1.5))))

integrand.one <- Vectorize(function(p) {g * product})
integrate(integrand.one, lower=0, upper=1)

integrand.two <- Vectorize(function(p) {g.pi.two * prod(p*dnorm(x, mean = 1, sd = 1) + (1-p)*dnorm(x, mean = 3, sd = sqrt(1.5)))})
integrate(integrand.two, lower=0, upper=1)

loglike <- function(p) { 
  alpha <- p[1]
  beta <- p[2]
  L <- log((o_ring$Failure*(exp(alpha+beta*o_ring$Temp))+(1-o_ring$Failure))/(1 + exp(alpha + beta*o_ring$Temp)))
  return(-sum(L))
}

estimates <- nlm(loglike, c(0,0))$estimate
estimates
