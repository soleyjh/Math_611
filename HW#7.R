# JAMES HYER SOLEY
# HOMEWOKR #7
# 10/29/2014

# QUESTION #2b
set.seed(100)
FUN.GMM <- function(){
  y <- runif(1)
  if (y <= .3){
    x <- rnorm(1, mean = 1, sd = 1)
  }
  else{
    x <- rnorm(1, mean = 3, sd = sqrt(1.5))
  }
  return(x)
}

samples <- sapply(1:100, function(x) FUN.GMM())
hist(samples)

# QUESTION #2c
set.seed(100)
likelihood <- function(p,x = samples) {
  L <- prod((p*dnorm(x, mean = 1, sd = 1) + (1-p)*dnorm(x, mean = 3, sd = sqrt(1.5))))
  return(L)
}

p <- seq(0,1,.001)
samples2 <- sapply(p, function(p) likelihood(p=p))
plot(p,samples2)

# QUESTION #2d
g.pi <- 1
x <- samples
likelihood.one <- expression(prod(p*dnorm(x,1,1)+(1-p)*dnorm(x,3,sqrt(1.5))))

integrand.one <- function(p){eval(likelihood.one)}
likelihood.prior.one <- sapply(p, function(x) integrand.one(x))
p.one <- integrate(Vectorize(integrand.one), lower=0, upper=1)$value
p.one

plot(p,likelihood.prior.one)

# QUESTION #2e(a)
x <- samples
g.two <- expression(6*p^5)
likelihood.two <- expression(prod(p*dnorm(x,1,1)+(1-p)*dnorm(x,3,sqrt(1.5))))

integrand.two <- function(p){eval(g.two)*eval(likelihood.two)}
likelihood.prior.two <- sapply(p, function(x) integrand.two(x))
p.two <- integrate(Vectorize(integrand.two), lower=0, upper=1)$value
p.two

plot(p,likelihood.prior.two)

#2e(b)
x <- samples
g.three <- expression(6*(1-p)^5)
likelihood.three <- expression(prod(p*dnorm(x,1,1)+(1-p)*dnorm(x,3,sqrt(1.5))))

integrand.three <- function(p){eval(g.three)*eval(likelihood.three)}
likelihood.prior.three <- sapply(p, function(x) integrand.three(x))
p.three <- integrate(Vectorize(integrand.three), lower=0, upper=1)$value
p.three

plot(p,likelihood.prior.three)


#QUESTION #3b
set.seed(100)
o_ring <- read.table("~/Math 611/HW#7/o_ring.txt", header=T, quote="\"")
likelihood <- function(p){
  alpha <- p[1] 
  beta <- p[2]
  L <- -sum(log((o_ring$Failure*(exp(alpha+beta*o_ring$Temp))+(1-o_ring$Failure))/(1+exp(alpha + beta*o_ring$Temp))))
  return(L)
}
  
mle <- nlm(likelihood, c(0,0))$estimate
mle

#QUESTION #3c
likelihood <- function(alpha,beta) {
  L <- prod((o_ring$Failure*(exp(alpha+beta*o_ring$Temp))+(1-o_ring$Failure))/(1 + exp(alpha + beta*o_ring$Temp)))
  return(L)
}

FUN.MH <- function(steps){
  
  # Initialize Alpha, Beta, Alpha_est, Beta_est
  alpha = 1
  beta = -1
  
  nu.i <- likelihood(alpha, beta)
  
  results <- data.frame("step"= 0, "alpha"= alpha, "beta"= beta, "likelihood" = nu.i, "measure" = alpha - beta)
    
  # Assume the Priors are Uniformly Distributed so they Cancel
  # g.alpha = pnorm(1, mean = 15.043, sd = 5)
  # g.beta = pnorm(1, mean = -.232, sd = 5)
    
  for(i in 1:steps){
    # Pick Alpha.Star and Beta.Star
    alpha.star <- alpha + rnorm(1)
    beta.star <- beta + rnorm(1)
    
    #Test Valid Propsal or not
    U <- runif(1)
    nu.j <- likelihood(alpha.star,beta.star)
    ratio <- exp(log(nu.j) - log(nu.i))
    
    #If valid make alpha.star = alpha and beta.star = beta
    if(U < min(ratio, 1) ) {
      alpha <- alpha.star
      beta <- beta.star
      nu.i <- nu.j
      results <- rbind(results, data.frame("step"=i, "alpha"=alpha, "beta"=beta, "likelihood"=nu.i, "measure"= alpha-beta))
      
    } 
    else {
      alpha <- alpha
      beta <- beta
      nu.i <- nu.i
      results <- rbind(results, data.frame("step"=i, "alpha"=alpha, "beta"=beta, "likelihood"=nu.i, "measure"= alpha-beta))

    }
    
  }
  return(results)
}

output <- FUN.MH(50000)
plot(output$measure)
plot(output$alpha)
plot(output$beta)

#Bayesian Estimates

failure <- function(alpha,beta, temp=65){
  return(exp(alpha+beta*temp)/(1+exp(alpha+beta*temp)))
}

fails.65 <- apply(output[10000:50000,c("alpha","beta")],1,function(x) failure(x['alpha'], x['beta']))
plot(density(fails.65))

failure <- function(alpha,beta, temp=45){
  return(exp(alpha+beta*temp)/(1+exp(alpha+beta*temp)))
}

fails.45 <- apply(output[10000:50000,c("alpha","beta")],1,function(x) failure(x['alpha'], x['beta']))
plot(density(fails.45))