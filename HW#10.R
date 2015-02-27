# JAMES HYER SOLEY
# HOMEWORK #10
# 11/16/2014

# QUESTION 2ai
# PULL MLE FUNCTION FROM HW#7
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

# Package numDeriv downloaded
library("numDeriv", lib.loc="~/R/win-library/3.1")

hessian <- hessian(func=likelihood,c(15.043,-.232))
hessian.inv <- solve(hessian)

alpha.ci.left <- 15.043 - 1.96*sqrt(hessian.inv[1,1])
alpha.ci.right <- 15.043 + 1.96*sqrt(hessian.inv[1,1])
alpha.ci <- cbind(alpha.ci.left,alpha.ci.right)
  
beta.ci.left <- -.232 - 1.96*sqrt(hessian.inv[2,2])
beta.ci.right <- -.232 + 1.96*sqrt(hessian.inv[2,2])
beta.ci <- cbind(beta.ci.left,beta.ci.right)

alpha.ci
beta.ci

# QUESTION 2aii
Ri <- o_ring$Failure
Ti <- o_ring$Temp
df.likelihood <- function(alpha, beta){
  df.alpha <- sum((Ri - 1)   / (Ri * exp(alpha + beta * Ti) - Ri +1) + (1  / (exp(alpha + beta*Ti)+1)))
return(cbind(df.alpha,df.beta))
  df.beta <-  sum((Ri*Ti-Ti) / (Ri* exp(alpha + beta * Ti) - Ri +1) + (Ti / (exp(alpha + beta*Ti)+1)))
}

x <- df.likelihood(15.043,-.232)
x.trans <- (t(x))

approx <- x.trans %*% x
approx.inv <- ginv(approx)

alpha.ci.left <- 15.043 - 1.96*sqrt(approx.inv[1,1])
alpha.ci.right <- 15.043 + 1.96*sqrt(approx.inv[1,1])
alpha.ci <- cbind(alpha.ci.left,alpha.ci.right)

beta.ci.left <- -.232 - 1.96*sqrt(approx.inv[2,2])
beta.ci.right <- -.232 + 1.96*sqrt(approx.inv[2,2])
beta.ci <- cbind(beta.ci.left,beta.ci.right)

alpha.ci
beta.ci

# QUESTION 2b
plot(density(Ti))
Ti.mean <- mean(Ti)
Ti.sd <- sd(Ti)

Ti.sim <- rnorm(100000,Ti.mean,Ti.sd)
P.Ri.sim <- exp(mle[1]+mle[2]*Ti.sim)/(1+exp(mle[1]+mle[2]*Ti.sim))
Ri.sim <- runif(100000)<P.Ri.sim


f.b <- function(X){
  X[1]-> R
  X[2]-> Ti
  (R*exp(mle[1]+mle[2]*Ti)+(1-R))/(1+exp(mle[1]+mle[2]*Ti))}


Z.sim <- sapply(1:100000, function(x) grad(f.b,c(Ri.sim[x],Ti.sim[x]))/f.b(c(Ri.sim[x],Ti.sim[x])))
ZZt.sim <- sapply(1:100000, function(x) Z.sim[,x]%*%t(Z.sim[,x]))
sd.ab.b <- sqrt(c(mean(ZZt.sim[1,]),mean(ZZt.sim[4,])))

alpha.ci.left <- 15.043 - 1.96*sd.ab.b[1]
alpha.ci.right <- 15.043 + 1.96*sd.ab.b[1]
alpha.ci <- cbind(alpha.ci.left,alpha.ci.right)

beta.ci.left <- -.232 - 1.96*sd.ab.b[2]
beta.ci.right <- -.232 + 1.96*sd.ab.b[2]
beta.ci <- cbind(beta.ci.left,beta.ci.right)

alpha.ci
beta.ci

# QUESTION 2c
likelihood.boot <- function(p){
  alpha <- p[1] 
  beta <- p[2]
  L <- -sum(log((bootstrap$Failure*(exp(alpha+beta*bootstrap$Temp))+(1-bootstrap$Failure))/(1+exp(alpha + beta*bootstrap$Temp))))
  return(L)
}

matrix <- matrix(0,nrow=10000,ncol=2)
for (i in 1:10000){
  bootstrap <- o_ring[sample(nrow(o_ring),nrow(o_ring)*10,replace=TRUE),]
  temp.mle <- nlm(likelihood.boot, c(0,0))$estimate
  se.temp <- temp.mle[1]-mle[1]
  se.failure <- temp.mle[2]-mle[2]
  error <- cbind(se.temp,se.failure)
  matrix[i,] <- cbind(error)  
  }

hist(matrix[,1], xlab="Alpha*", main = "Histogram of Alpha* - Alpha_Mle")
hist(matrix[,2], xlab="Beta", main = "Histogram of Beta* - Beta_Mle")


# QUESTION #3a
x.sim <- function(rolls, alpha) {
  #we say F = 0 and C = 1
  x <- 0
  path <- x
  for (i in 1:rolls) {
    U <- runif(1)
    x.new <- ifelse(U <= alpha, x <- abs(x-1), x <- x)
    path <- c(path,x.new)
  }
  
  return(path)
}

y.sim <- function(x.state) {
  ifelse(x.state == 0, y <- sample(1:6,1), y <- sample(c(6, rep(1:5,3)),1))
  return(y)
}

chain.sim <- function(rolls, alpha) {
  x<-x.sim(rolls, alpha)
  chain <- sapply(x, function(x) y.sim(x))
  chain.final <- data.frame("X.State" = x, "Y.Roll" = chain)
  return(chain.final)
}

#M-H algorithm
G <- function(state, roll) {
  if(state == "F") {
    return(1/6)
  } else if ( roll == 6) {
    return(1/16)
  } else {
    return(3/16)
  }
}

P.change <- function(alpha, s.curr, s.prev) {
  return(ifelse(s.curr == s.prev, 1-alpha, alpha))
}

CasinoMCMC <- function(time) {
  current <- chain.sim(1000, .02)
  path <- as.data.frame(cbind("Time" = 0, current))
  i <- 1
    
  while ( i <= time) {
    # Make a Proposal
    proposal <- chain.sim(1000, .02)
    
    # Check Validity of Proposal
    nu.i <- sum(log(mapply(G, current$X.State[1:1000], current$Y.Roll[1:1000]))) + sum(log(mapply(P.change, .02, current$X.State[1:1000], current$X.State[2:1001]))) + log(G(current$X.State[1001], current$Y.Roll[1001]))
    nu.j <- sum(log(mapply(G, proposal$X.State[1:1000], proposal$Y.Roll[1:1000]))) + sum(log(mapply(P.change, .02, proposal$X.State[1:1000], proposal$X.State[2:1001]))) + log(G(proposal$X.State[1001], proposal$Y.Roll[1001]))
    
    # Generate Random Uniform Variable for Valid Accept / Reject
    U <- runif(1)
    
    # If Valid Proposal nu.i = nu.j & current = proposal
    # Else create a new proposal
    
    if (U < min(1,exp(nu.i-nu.j))) {
      current <- proposal
      nu.i <- nu.j
      path <- rbind(path, as.data.frame(cbind("Time" = i, current)))
    } else {
      path <- rbind(path, as.data.frame(cbind("Time" = i, current)))
    }
    i <- i + 1
  }
    path$row_num <- 1:nrow(path)
    path$prob_cheat <- cumsum(path$X.State)/path$row_num
  return(path$prob_cheat)
}



results <- mclapply(1:100100, function(i) CasinoMCMC(100) , mc.cores = 4 )

library(Parallel)
cl <- CasinoMCMC(100)
registerDoParallel(cl)