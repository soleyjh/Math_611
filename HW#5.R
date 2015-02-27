#James Hyer Soley
#STOCHASTIC SIMULATION MATH 611
#10/14/14

#QUESTION #2b
#set starting parameters
hiv <- read.table("~/Math 611/HW#5/HIV.txt", header=T, quote="\"")
params <- data.frame(alpha = .2, beta =.5, mu = 5, lambda = 10)

loglike <- function(n = hiv$frequency, i = hiv$encounters, params = data.frame(alpha = .01, beta =.02, mu = .7, lambda = .3)) {
  #using dpois() and 0! is 1, the factorials are integrated into the loglikelihood function for better readibility and to reduce human error
  pi <- ifelse(i == 0, (params$alpha + params$beta*dpois(i, params$mu) + (1-params$alpha-params$beta)*dpois(i,params$lambda)),(params$beta*dpois(i, params$mu) + (1-params$alpha-params$beta)*dpois(i,params$lambda)))
  fxn <- sum(n * log(pi))
  return(fxn)
}

EM.hiv <- function(n = hiv$frequency, i=hiv$encounters, params = data.frame(alpha = .01, beta =.02, mu = .7, lambda = .3)) {
  iter <- 2
  L.theta <- loglike(params)
  path <- rbind(list(iter = 0, alpha = 0, beta = 0, mu = 0, lambda = 0, L.theta = 0),list(iter=0, alpha = params$alpha, beta = params$beta, mu = params$mu, lambda = params$lambda, L.theta = L.theta))
  
  while (abs(path[iter,]$L.theta - path[iter-1,]$L.theta) > 0) {
    
    #E step
    z <- params$alpha/(params$alpha + params$beta*dpois(0, params$mu) + (1-params$alpha-params$beta)*dpois(0,params$lambda))
    t <- ifelse(i == 0, ((params$beta*dpois(i, params$mu)) /(params$alpha + params$beta*dpois(i, params$mu) + (1-params$alpha-params$beta)*dpois(i,params$lambda))),((params$beta*dpois(i, params$mu))/(params$beta*dpois(i, params$mu) + (1-params$alpha-params$beta)*dpois(i,params$lambda))))
    p <- ifelse(i == 0, (((1 - params$alpha - params$beta)*dpois(i, params$lambda)) /(params$alpha + params$beta*dpois(i, params$mu) + (1-params$alpha-params$beta)*dpois(i,params$lambda))),(((1 - params$alpha - params$beta)*dpois(i, params$lambda))/(params$beta*dpois(i, params$mu) + (1-params$alpha-params$beta)*dpois(i,params$lambda))))
    
    
    #M step
    params$alpha <- n[1]*z/sum(n)
    params$beta <- sum(n*t/sum(n))
    params$mu <- sum(i*n*t)/sum(n*t)
    params$lambda <- sum(i*n*p)/sum(n*p)
    
    L.theta <- loglike(params)
    iter <- iter + 1
    path <- rbind(path,list(iter = (iter-2), alpha = params$alpha, beta = params$beta, mu = params$mu, lambda = params$lambda, L.theta = L.theta))
  }
  
  return(as.list(path[2:nrow(path),]))
}

path <- EM.hiv(params)

library(ggplot2)
qplot(x=as.numeric(path[,1]), y=as.numeric(path[,6]), xlab = "Iteration", ylab="L.theta")


#QUESTION #3a
set.seed(100)
FUN.MARKOV= function(Time, ro){
  x=c()
  P=matrix(c(1-3*ro,ro,2*ro,ro,1-5*ro,4*ro,2*ro,4*ro,1-6*ro),nrow=3)
  x[1]=1
  for (i in 1: Time){
    u=runif(1,0,1)
    temp= u-cumsum(P[x[i],]) 
    x[i+1]= which(temp<=0)[1]
  }
  return(x)
}
FUN.MARKOV(6,0.3)

#QUESTION #2bi
FUN.MARKOV.1= function(ro){
  P=matrix(c(1-3*ro,ro,2*ro,ro,1-5*ro,4*ro,2*ro,4*ro,1-6*ro),nrow=3)
  pi=eigen(P)$vectors[,which(eigen(P)$values==max(eigen(P)$values))]
  return(pi)
}

FUN.MARKOV.1(0.1)
FUN.MARKOV.1(0.01)
FUN.MARKOV.1(0.0001)

#QUESTION #2bii
FUN.MARKOV.2=function(Time,ro){
  pi=c()
  FUN.TIME = FUN.MARKOV(Time,ro)
  pi[1]=sum(FUN.TIME==1)/Time
  pi[2]=sum(FUN.TIME==2)/Time
  pi[3]=sum(FUN.TIME==3)/Time
  return(pi)
  
}
FUN.MARKOV.2(10000,0.1)
FUN.MARKOV.2(10000,0.01)
FUN.MARKOV.2(10000,0.001)

#QUESTION #2biii
FUN.MARKOV.3=function(Time,ro){
  P=matrix(c(1-3*ro,ro,2*ro,ro,1-5*ro,4*ro,2*ro,4*ro,1-6*ro),nrow=3)
  Pt=P
  for (i in 1:Time){
    Pt=Pt%*%P
  }
  pi=Pt[,1]
  return(pi)
  
}

FUN.MARKOV.3(10000,0.1)
FUN.MARKOV.3(10000,.01)
FUN.MARKOV.3(10000,0.0001)