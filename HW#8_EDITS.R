#Elaine Code


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

sim.1000 <- chain.sim(1000,.02)
plot(x = 50000:200001, y = sim.200k$X.State[50000:200001], type = "l")

(1/150000) * sum(sim.200k$X.State[50000:200000])

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

mh.sampler <- function(time) {
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
  return(path)
}

# Run 10,000 Times is about 1 hour
x <- mh.sampler(10000)

# Drop the burn time (roughly 1/2 of samples, based on class discussions)
non.burn <- sum(x$X.State[5005501:10011001]) 
estimate <- non.burn/(10011001-5005501)
estimate