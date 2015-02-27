# JAMES HYER SOLEY
# HOMEWORk #8
# 11/05/2014

# QUESTION #2a
set.seed(100)
samples <- function(n,rate){
  x <- cumsum(rexp(n,rate))
  return(x)
}
samples(45,1)

# QUESTION #2b
# Create Function Lambda
lambda <- function(t){
  rate <- 5e-04 + 7.5858e-05*exp(log(1.09144)*t)
}

# Create Death Function
death <- function(t){
  x <- samples(10000, lambda(t))
  max <- lambda(t)
  proposals <- x[which(x < t)]
  alive <- TRUE

  for (i in 1:length(proposals)){
    u <- runif(1)    
    if (alive == TRUE && u < lambda(proposals[i])/max){
      alive == FALSE
      return(proposals[i])
    }
    else {
      i = i + 1
    }
  }  
}

death(100)
death(100)
death(100)
death(100)
death(100)
death(100)
death(100)
death(100)
death(100)


# QUESTION #2d 
# Integrate the function for k

integrand <- function(x) {return((5*10^(-4)+7.5858*10^(-5)*exp(log(1.09144)*x)))}
exponent <- integrate(Vectorize(integrand), 60, 90)
P.90 <- integrate(Vectorize(integrand), 0, 90)
P.60 <- integrate(Vectorize(integrand), 0, 60)

# P(T>90 | T >60) = P(T>90) = [1 - P(T<90)] / [1 - P(T<60)]
integrate.death <- exp(-P.90$value)/exp(-P.60$value)
integrate.death

# Monte Carlo
deaths <- sapply(1:100000, function(x) death(90))
gt.90 <- sum(length(deaths[which(deaths > 90)]))
gt.60 <- sum(length(deaths[which(deaths > 60)]))

mc.death <- gt.90/gt.60
mc.death

# QUESTION #3a
x.t <- function(rolls, alpha) {
  x <- 1
  path <- x
  for (i in 1:rolls) {
    U <- runif(1)
    x.new <- ifelse(U <= alpha, x <- -x, x <- x)
    path <- c(path,x.new)
  }
    final <- ifelse(path == 1, x <- "F", x <- "C")
  return(final)
}

y.t <- function(x.state) {
  ifelse(x.state == "F", y <- sample(1:6,1), y <- sample(c(6, rep(1:5,3)),1))
  return(y)
}

chain.sim <- function(rolls, alpha) {
  x<-x.t(rolls, alpha)
  chain <- sapply(x, function(x) y.t(x))
  chain.final <- data.frame("X.State" = x, "Y.Roll" = chain)
  return(chain.final)
}

# QUESTION #3b
sim.1000 <- chain.sim(1000, .02)
head(sim.1000)

# QUESTION #3d
#M-H algorithm
G <- function(state, roll) {
  if(state == "F") {
    return(1/6)
  } else if (roll == 6) {
    return(1/16)
  } else {
    return(3/16)
  }
}

P.change <- function(alpha, s.curr, s.prev) {
  return(ifelse(s.curr == s.prev, 1-alpha, alpha))
}

mh.sampler <- function(time) {
  current <- chain.sim(1000, .02)
  path <- as.data.frame(cbind("Time" = 0, current))
  i <- 1
  
  while ( i <= time) {
    proposal <- chain.sim(1000, .02)
    
    nu.i <- sum(log(mapply(G, current$X.State[1:1000], current$Y.Roll[1:1000]))) + sum(log(mapply(P.change, .02, current$X.State[1:1000], current$X.State[2:1001]))) + log(G(current$X.State[1001], current$Y.Roll[1001]))
    nu.j <- sum(log(mapply(G, proposal$X.State[1:1000], proposal$Y.Roll[1:1000]))) + sum(log(mapply(P.change, .02, proposal$X.State[1:1000], proposal$X.State[2:1001]))) + log(G(proposal$X.State[1001], proposal$Y.Roll[1001]))
    
    U <- runif(1)
    
    if (U < min(1,exp(nu.i-nu.j))) {
      current <- proposal
      nu.i <- nu.j
      path <- rbind(path, as.data.frame(cbind("Time" = i, current)))
    } else {
      path <- rbind(path, as.data.frame(cbind("Time" = i, current)))
    }
    i <- i + 1
    print(i)
  }
  return(path)
}

