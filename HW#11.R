# JAMES HYER SOLEY
# HOMEWORK #11
# 12/6/2014

# QUESTION 1a
samp.brown <- function(t){
  x.t <- cumsum(rnorm(t,sqrt(.01)))
  return (x.t)
}

x.1 <- samp.brown(10)
plot(x.1, main = "Path 1")

x.2 <- samp.brown(10)
plot(x.2, main = "Path 2")

x.3 <- samp.brown(10)
plot(x.3, main = "Path 3")

#QUESTION 1ci
random.walk <- function(n){
  x.t <- cumsum(sample(c(-1,1),n, replace = TRUE))
  return(x.t)
}

plot(random.walk(1000000), xlim=c(0,1000000), ylim=c(-1000000,1000000))
plot(random.walk(1000000), xlim=c(0,1000000), ylim=c(-3*10^3,+3*10^3))


#QUESTION 2a
election.data <- read.delim("~/Math 611/HW#11/2012_Election_Polling_Data.txt")

n <- nrow(election.data)
elec.ci <- c(-1,-1, 0,-1)
elec.ui <- cbind(c(1,-1,0,0), c(0,0,1,-1))

election.model <- function(mu_sigma) {
  actual.Obama <- 0.514
  actual.Romney <- 1-actual.Obama
  y.i <- rep(NA, n)
  
  for (x in 1:n) {
    normalized <- sqrt(actual.Obama * actual.Romney/election.data[x,2])
    y.i[x] <- rnorm(1, actual.Obama, normalized) + rnorm(1, mu_sigma[1], mu_sigma[2])
  }
  return(-sum(log(y.i)))
}

set.seed(611)
constrOptim(c(-0.1, 0.1), election.model, ui=elec.ui, ci=elec.ci, grad=NULL)

#QUESTION 2b
samples <- 1e3
mle.samples <- matrix(NA, samples, 2)
colnames(mle.samples) <- c("mu", "sigma")
elec.ci.boot <- c(-0.5,-0.5,0,-1)
elec.ui.boot <- cbind(c(1,-1,0,0), c(0,0,1,-1))

for (y in 1:samples) {
  X <- sample(1:n, n, replace=TRUE)
  election.sample <- data.frame(election.data[X,])
  mle <- constrOptim(c(0,0.1), election.model, ui=elec.ui.boot, ci=elec.ci.boot, grad=NULL)
  mle.samples[y,1] <- mle$par[1]
  mle.samples[y,2] <- mle$par[2]
}

par(mfrow=c(2,2))
hist(mle.samples[,1], breaks=100, main="", xlab=paste(expression(mu)))
hist(mle.samples[,2], breaks=100, main="", xlab=paste(expression(sigma)))
plot(sort(mle.samples[,1]), ylab="Samples")
plot(sort(mle.samples[,2]), ylab="Samples")

mu.lower <- sort(mle.samples[,1])[0.025*samples]
mu.upper <- sort(mle.samples[,1])[0.975*samples]
sigma.lower <- sort(mle.samples[,2])[0.025*samples]
sigma.upper <- sort(mle.samples[,2])[0.975*samples]

#QUESTION 3
maturityDateT <- 30;
e <- 0.1;
std <- 0.30;
dT <- (e)^2/(std * 5)^2;

AsianOption <- function(t, delT = dT, s0 = 8, mu = -0.045, sdev = 0.30){
  
  runs <- t/delT + 1;
  strike <- 10;
  
  B <- cumsum(rnorm(n = runs, mean = 0, sd = sqrt(delT)));
  time <- seq(from = 0, to = 30, by = delT);
  Price <- function(i) s0 * exp(mu * time[i] + sdev * B[i]) - strike;
  out <- sapply(1:(runs-1), Price);
  
  over.strike <- function(i) max(0,out[i]);
  St <- sapply(1:(runs-1),over.strike);
  
  return(sum(dT * St)/maturityDateT);
}

C <- rep(NA, 2e4);
for (j in 1:length(C)){
  C[j] <- AsianOption(maturityDateT);
  cat(j, "\n"); 
}

hist(C, breaks = 5000, xlim = c(0,5), main = "Histogram of Price of Call", xlab = "Dollar Val");

avg.price <- sum(C)/j;
avg.price
