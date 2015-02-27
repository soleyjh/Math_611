# JAMES HYER SOLEY
# HOMEWORk #9
# 11/14/2014

# QUESTION #2a
# Read in Data
gdp <- read.table("~/Math 611/HW#9/gdp.txt", header=TRUE, quote="\"")

# Page 372 -- Parameters Given
alpha.one = 1.522
alpha.zero = -0.3577
p = .9049
q = .7550
sigma = 0.7690
phi.one = .014 
phi.two = -.058 
phi.three = -.247
phi.four=-.213

# From GDP Data Given
gpd <- gdp[,2]
gpd.log = 100*(log(gpd[2:267]) - log(gpd[1:266]))

# Define the Error Function 
error = function(z.vect, phi.one, phi.two, phi.three, phi.four, alpha.zero, alpha.one, sigma){
  dnorm((z.vect[5:266] - phi.one*z.vect[4:265]- phi.two*z.vect[3:264] -  
           phi.three*z.vect[2:263] - phi.four*z.vect[1:262]),0,sigma) 
}

# Define Transition Function
transition = function(states, p, q){
  c(sapply(1:(length(states)-1), function(i) ifelse(states[i]==1 & states[i+1]==1, p , ifelse(states[i]==0
                              & states[i+1]==0, q, ifelse(states[i] == 1 & states[i+1] == 0, 1-p,1-q)))))
}

# Define the NU Function
nu = function(trans.prob, error.prob){
  sum(log(c(trans.prob,error.prob)))
}

# MCMC Function
MCMC = function(step,gpd,states,phi.one,phi.two,phi.three,phi.four,alpha.zero,alpha.one,sigma){
  i = 1
  states = sample(c(0,1), 266, replace = TRUE)
  store.metr <- rep(0, step)
  store.states <- matrix(c(0), nrow = 266, ncol = step)
    
  while(i < step){
    flip <- sample(1:266,1, replace = FALSE)
    z.vect.old = -alpha.one * states - alpha.zero + gpd.log
    error <- error(z.vect.old, phi.one, phi.two, phi.three, phi.four, alpha.zero, alpha.one,sigma)
    trans <- transition(states, p, q)
    
    states.new <- states
    states.new[flip] <- -states.new[flip] + 1
    z.vect.new = -alpha.one * states.new - alpha.zero + gpd.log
    error.new <- error( z.vect.new, phi.one, phi.two, phi.three, phi.four, alpha.zero, alpha.one,sigma)
    trans.new <- transition(states.new, p, q)
    
    if(runif(1) < min(1, exp(nu(trans.new, error.new) - nu(trans, error)))){
      states <- states.new
      store.states[,i] <- states
      store.metr[i] <- sum(states)
      i = i + 1
    }
  }
  return(list("metric" =store.metr[-step], "states" = store.states, "final.state" = states))
}

# Initialize the states

new.run <- MCMC(10000,gpd,states,phi.one, phi.two, phi.three, phi.four, alpha.zero, alpha.one, sigma)
plot(new.run$metric, main = "Steady State")
plot(1-rowSums(new.run$states[,5000:10000])/5000, type = "l", main = "Recession Probability", xlab = "Time", ylab = "Probability @ S(t)" )
plot( (1-rowSums(new.run$states[,5000:10000])/5000) > .5 , type = "l", main = "Recessions > .5 Probability", xlab = "Time", ylab = "Probability @ S(t)")


# Queston #2b
# Update MCMC Function from above
MCMC.BAYESIAN = function(step,gpd){
  # Create the Holders and Variables needed to record Analysis
  i = 1
  store.metr <- rep(0, step)
  store.states <- matrix(c(0), nrow = 266, ncol = step)
  
  # Initialize All Parameters Arbitrarily
  alpha.zero <- 0
  alpha.one <- 1
  phi.one <- 0
  phi.two <- 0
  phi.three <- 0
  phi.four <- 0
  sigma <- 1
  p <- .5
  q <- .5
  
  # Calculate nu.i
  # Assume Uniform Priors so they cancel
  
  states.new <- states
  states.new[flip] <- -states.new[flip] + 1
  states = sample(c(0,1), 266, replace = TRUE)
  flip <- sample(1:266,1, replace = FALSE)
  z.vect = -alpha.one * states - alpha.zero + gpd.log
  error <- error(z.vect, phi.one, phi.two, phi.three, phi.four, alpha.zero, alpha.one,sigma)
  trans <- transition(states, p, q)

  #Form nu.i ratio
  nu.i <- exp(nu(trans.new, error.new) - nu(trans, error)
  
  while(i < step){
    alpha.zero.star <- alpha.zero + rnorm(1,0,.01)
    alpha.one.star <- alpha.one + rnorm(1,0,.01)
    phi.one.star <- phi.one + rnorm(1,0,.01)
    phi.two.star <- phi.two + rnorm(1,0,.01)
    phi.three.star <- phi.three + rnorm(1,0,.01)
    phi.four.star <- phi.four + rnorm(1,0,.01)
    sigma.star <- sigma + rnorm(1,0,.01)
    p.star <- p + rnorm(1,0,.01)
    q.star <- q + rnorm(1,0,.01)
    
    states.new <- states
    states.new[flip] <- -states.new[flip] + 1
    z.vect.new = -alpha.one.star * states.new - alpha.zero.star + gpd.log
    error.new <- error( z.vect.new, phi.one.star, phi.two.star, phi.three.star, phi.four.star, alpha.zero.star, alpha.one.star,sigma.star)
    trans.new <- transition(states.new, p.star, q.star)
    
    # Form nu.i and nu.j for Ratio
    nu.j <- exp(nu(trans.new, error.new) - nu(trans, error)
    ratio <- nu.i/nu.j
                
    if(runif(1) < min(1, ratio){
      states <- states.new
      store.states[,i] <- states
      store.metr[i] <- sum(states)
      nu.i <- nu.j
      
      # Change the 
      alpha.zero <- alpha.zero.star
      alpha.one <- alpha.one.star
      phi.one <- phi.one.star
      phi.two <- phi.two.star
      phi.three <- phi.three.star
      phi.four <- phi.four.star
      sigma <- sigma.four.star
      p <- p.star
      q <- q.star
      
      i = i + 1      
    }
    else 
      states <- states.new
      store.states[,i] <- states
      store.metr[i] <- sum(states)
      
  }
  return(list("metric" =store.metr[-step], "states" = store.states, "final.state" = states))
}


