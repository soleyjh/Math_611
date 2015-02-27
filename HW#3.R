#JAMES HYER SOLEY
#STOCHASTIC SIMULATION MATH 611
#9/24/2014


#QUESTION #2bi
set.seed(100)
FUN.GMM <- function(){
  y <- runif(1)
  if (y <= .7){
    x <- rnorm(1, mean = 0, sd = 1)
    }
  else{
    x <- rnorm(1, mean = 3, sd = 2)
    }
  return(x)
}

ans <- NULL
for (i in 1:100){
  ans[i] <- FUN.GMM()
}

hist(ans)

# QUESTION #2bii (easy function to start with)
# A steepest accent algorithm to find the maximum

M <- c(.2,.1,.1,.5,.3)
FUN.GRAD <- function(M) {
  Q <- M[1]*((2*pi*M[4])^-.5)*exp(-(X-M[2])^2/(2*M[4])) + (1-M[1])*((2*pi*M[5])^-.5)*exp(-(X-M[3])^2/(2*M[5]))
  p1 <- (sum((((2*pi*M[4])^-.5)*exp(-(X-M[2])^2/(2*M[4])) - ((2*pi*M[5])^-.5)*exp(-(X-M[3])^2/(2*M[5])))/Q)
  mu1 <- sum((M[1]*((2*pi*M[4])^-.5)*exp(-(X-M[2])^2/(2*M[4]))*2*(X-M[2]))/Q)
  mu2 <- sum((1-M[1])*((2*pi*M[5])^-.5)*exp(-(X-M[3])^2/(2*M[5]))*2*(X-M[3])/Q)
  sig1 <- sum((M[1]/(sqrt(2*pi)*Q*M[4]^4))*exp(-(X-M[2])^2/(2*M[4]))*(M[4]^2.5*(-X+M[2])-.5*M[4]^2.5 +.5*M[4]^1.5*(X-M[2])^2))
  sig2 <- sum(((1-M[1])/(sqrt(2*pi)*Q*M[5]^4))*exp(-(X-M[3])^2/(2*M[5]))*(M[5]^2.5*(-X+M[3])-.5*M[5]^2.5 +.5*M[5]^1.5*(X-M[3])^2))
  grad <- c(p1, mu1, mu2, sig1, sig2)
  return(grad) 
}

# Calculate the Gradient norm
FUN.GRAD_NORM <- function(grad) {
  gradnorm = sqrt(sum(grad^2))
  return(gradnorm)
}

FUN.STEEP <- function(M, step){
  while(gradnorm > .0001){
    slope <- grad/gradnorm
    M <- M+slope*step
    FUN.GRAD(M)
    FUN.GRAD_NORM(gradnorm)
  }
}

step <- .001
FUN.STEEP(M, step)
