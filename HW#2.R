#JAMES HYER SOLEY
#STOCHASTIC SIMULATION MATH 611
#9/17/2014

#QUESTION # 3a
FUN.QNORM <- function(p,start,step){
  norm <- function(x){1/sqrt(2*pi)*exp(-x^2/2)}
  area <- 0
  x <- start + step
  while (area < p){
    integrand <- integrate(norm,start,x)
    area <- integrand$value
    x = x + step
    }
    return(x - step)
}

#QUESTION # 3a FIRST PART
p <- runif(1)
system.time(FUN.QNORM(p,-100,.001))
system.time(qnorm(p))

#QUESTION # 3b FIRST PART
p <- runif(1)
FUN.QNORM(p,-100,.001)
qnorm(p)

#QUESTION #3b SECOND PART
set.seed(100)
ans <- NULL
for (i in 1:10000){
  p <- runif(1)
  ans[i] <- FUN.QNORM(p,start = -5,step = .001)
  if(i==5000){
  print('Observation 5000:')
  print(qnorm(p))
  print(ans[i])
  }
}

hist(ans)

#QUESTION #4
#Initialize Both A and P for Our Model

set.seed(100)
a <- 1:1000000
p <- .00001*(.99999)^(a-1)
sum(p)

#Create Function to Sample Discretely
FUN.RGEOM <- function(p){
  x <- runif(1)
  prob <- 0
  for (i in 1:length(p)){
    if(prob > x)
      break
    prob <- sum(p[1:i])
  }
return(i)
}

# Sample 10 Times
ans <- NULL
system.time(for (j in 1:10){
  ans[j] <- FUN.RGEOM(p)
})

ans
r_ans <- rgeom(10,.00001)

#Question 3c
#C is 2 determined by some simple math and knowledge of distributions
#N is the desired number of samples
genNorm <- function(N){
  nvars <- c()
  
  while(length(nvars) < N){
    u <- 2*runif(1); 
    expvars <- runif(1)
    
    if (expvars > .5) {
      y = rexp(1,rate=1)
    } else {
      y = -rexp(1, rate=1)
    }
    
    if(u < exp((y^2/2)-y)/sqrt(2*pi)) {
      nvars <- c(nvars,y)
    }
    
  }
  return(nvars)
}

vars <- genNorm(10000)
hist(vars)
qqnorm(vars)

rnormals <- rnorm(10000)
hist(rnormals)
qqnorm(rnormals)