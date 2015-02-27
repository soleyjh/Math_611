#JAMES HYER SOLEY
#STOCHASTIC SIMULATION MATH 611
#10/04/2014

#QUESTION #2a
set.seed(100)

S <- matrix(c(2,1,1,10),2,2)
Q <- eigen(S)

sampler <- function(Q, n) {
  x <- y <- matrix(NA,2,n)
  
  for (i in 1:n) {
    x[,i] <- Q$vectors %*% rnorm(2,1,sqrt(Q$values))
  }
  return(x)
}

test.x <- sampler(Q,1000000)
plot(test.x[1,],test.x[2,], xlim = c(-15, 15),ylim=c(-20,20), pch=16, cex=0.5)

#QUESTION #3aii
set.seed(100)

hope <- read.table("~/Math 611/HW#4/hope.txt", header=T, quote="\"")

heights <- hope$Height
means <- c(49,90)

norm <- function(x,mean) {
  norm <- sqrt(sum((x-mean)^2))
  return(norm)
}

k.means <- function(testdata, means) {
  data <- cbind(testdata, c(1:length(testdata)))
  test.init <- 1
  test.final <- 0
  mu1 <- means[1]
  mu2 <- means[2]
  means <- rbind(c(0,0),c(mu1,mu2))
  distances <- NULL
  i <- 2 
  
  while(abs(means[i,1] - means[i-1,1]) > .01 | abs(means[i,2] - means[i-1,2]) > .01) {
    norm.mus <- sum(sapply(data[,1], norm, mean = data[,2])^2)
    distances <- c(distances,norm.mus)
    
    #assign x's based on mus
    data.1 <- sapply(data[,1], norm, mean = mu1 )
    data.2 <- sapply(data[,1], norm, mean = mu2 )
    data[,2] <- ifelse(data.1^2 < data.2^2, mu1, mu2)
    
    #calculate new means based on classes
    mu1 <- mean(subset(data[,1], data[,2] == mu1 ))
    mu2 <- mean(subset(data[,1], data[,2] == mu2 ))
    means <- rbind(means,c(mu1,mu2))
    
    i <- i+1
  }
  distances <- c(distances,norm.mus)
  clusters <- ifelse(data[,2] == mu1, 1, 2)
  results <- list(means = c(mu1, mu2), clusters = as.vector(clusters), meanchanges = means, distances = distances)
  return(results)
}

mymeans <- k.means(heights, means)
rmeans <- kmeans(heights, centers=2)
cbind(mymeans$cluster,rmeans$cluster)
plot(x=mymeans$distances, type ="l", ylab = "Value of f(x)", xlab = "Iteration")


#QUESTION #3aiii
diffs <- cbind(ifelse(mymeans$cluster == hope$Gender, 0, 1),ifelse(rmeans$cluster == hope$Gender, 1, 0))
1 - sum(diffs[,1])/100
1 - sum(diffs[,2])/100


#QUESTION 3C
#EM Algorithm
set.seed(100)
hope <- read.table("~/Math 611/HW#4/hope.txt", header=T, quote="\"")

heights <- hope$Height
loglike <- function(x,mu.one,mu.two,sigma.one,sigma.two,p) {
  fxn <- sum( log(p*(1/sqrt(2*pi*sigma.one))*exp(-(x-mu.one)^2/(2*sigma.one)) + (1-p)*(1/sqrt(2*pi*sigma.two))*exp(-(x-mu.two)^2/(2*sigma.two))) )
  return(fxn)
}

#EM
EM.est <- function() {
  #starting variance = variance of the entire sample 
  #starting mu = two random values from the sample
  x <-heights
  mu1 <- sample(heights, 1, replace=FALSE)
  mu2 <- sample(heights, 1, replace=FALSE)
  sig1 <- 16
  sig2 <- 17
  p <- .5
  log <- loglike(x, mu1, mu2, sig1, sig2, p)
  path <- rbind(c(0,0,0,0,0,0),c(mu1, mu2, sig1, sig2, p, log))
  i <- 2
  
  #sufficiently small = .01
  while (abs(path[i,1] - path[i-1,1]) > .01 | abs(path[i,2] - path[i-1,2]) > .01 ) {
    
    #expectation E
    r1 <- p*(sqrt(2*pi*sig1)^(-1))*exp((-(x-mu1)^2)*(2*sig1)^(-1)) / (p*(sqrt(2*pi*sig1)^(-1))*exp((-(x-mu1)^2)*(2*sig1)^(-1)) + (1-p)*(sqrt(2*pi*sig2)^(-1))*exp((-(x-mu2)^2)*(2*sig2)^(-1)))
    r2 <- (1-p)*(sqrt(2*pi*sig2)^(-1))*exp((-(x-mu2)^2)*(2*sig2)^(-1)) / (p*(sqrt(2*pi*sig1)^(-1))*exp((-(x-mu1)^2)*(2*sig1)^(-1)) + (1-p)*(sqrt(2*pi*sig2)^(-1))*exp((-(x-mu2)^2)*(2*sig2)^(-1)))
    
    #maximization M
    mu1 <- sum((r1) * x)/sum((r1))
    mu2 <- sum(r2 * x)/sum(r2)
    sig1 <- sum((r1) * (x-mu1)^2)/sum(r1)
    sig2 <- sum(r2 * (x-mu2)^2)/sum(r2)
    p <- 1/100 * sum(r2)
    p2 <- 1/100 * sum(r1)
    
    log <- loglike(x, mu1, mu2, sig1, sig2, p)
    i <- i + 1
    path <- rbind(path,c(mu1, mu2, sig1, sig2, p, log))
  }
  
  return(path)
}

#QUESTION 3ci
#Graph of log likelihood
path <- EM.est()
plot(path[2:nrow(path),6], ylim=c(-285,-280), xlab="Iteration", ylab="Loglikelihood")

#QUESTION 3cii
x <- hope$Height
mu1 <- path[11,1]
mu2 <- path[11,2]
sig1 <- path[11,3]
sig2 <- path[11,4]
p <- path[11,5]

r1 <- p*(sqrt(2*pi*sig1)^(-1))*exp((-(x-mu1)^2)*(2*sig1)^(-1)) / (p*(sqrt(2*pi*sig1)^(-1))*exp((-(x-mu1)^2)*(2*sig1)^(-1)) + (1-p)*(sqrt(2*pi*sig2)^(-1))*exp((-(x-mu2)^2)*(2*sig2)^(-1)))
r2 <- (1-p)*(sqrt(2*pi*sig2)^(-1))*exp((-(x-mu2)^2)*(2*sig2)^(-1)) / (p*(sqrt(2*pi*sig1)^(-1))*exp((-(x-mu1)^2)*(2*sig1)^(-1)) + (1-p)*(sqrt(2*pi*sig2)^(-1))*exp((-(x-mu2)^2)*(2*sig2)^(-1)))
class <- ifelse (r1 > r2, 1, 2)
diffs <- sum(ifelse(class == hope$Gender, 1, 0))
1-sum(diffs)/100