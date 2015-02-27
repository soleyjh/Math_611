#James Hyer Soley
#STOCHASTIC SIMULATION MATH 611
#10/19/14

#QUESTION 2a
#Initialize the 100 x 100 Matrix

FUN.COIN <- function(step){
  graph <- matrix(0,100,100)
  for (i in 1:step){
    x <- trunc(runif(1,1,101))
    y <- trunc(runif(1,1,101))
    flip <- trunc(runif(1,0,2))
    
    if(flip){
      change <- TRUE
      if(x+1 <= 100){
        if(graph[x+1,y] == 1){change<-FALSE}
      }
      if(x >= 2 & change){
        if(graph[x-1,y] == 1){change<-FALSE}
      }
      if(y+1 <= 100 & change){
        if(graph[x,y+1] == 1){change<-FALSE}
      }        
      if(y >= 2 & change){
        if(graph[x,y-1] == 1){change<-FALSE}
      }
      if(change){graph[x,y]<-1}
    }
    else graph[x,y]<-0
  }
return(sum(graph))
}

#Plot X by time
x <- FUN.COIN(40000)
plot(x)

#Generate 1000 Samples of X
samp <- NULL
for (i in 1:1000){
  samp[i]<-(FUN.COIN(20000)/10000)
}
hist(samp)


#Generate 1000 samples of x and normalize
fx <- NULL
FUN.COIN1 <- function(step){
  #Initialize for #MH Ratio Test Proposal
  graph <- matrix(0,100,100)
  nu.i < - 1
  nu.j < - sum(graph)^2
  
  for (i in 1:step){
    
    u <- runif(1)
    if (u < min(1,nu.j/nu.i)){
    
    x <- trunc(runif(1,1,101))
    y <- trunc(runif(1,1,101))
    flip <- trunc(runif(1,0,2))
    
  if(flip){
      change <- TRUE
      if(x+1 <= 100){
        if(graph[x+1,y] == 1){change<-FALSE}
      }
      if(x >= 2 & change){
        if(graph[x-1,y] == 1){change<-FALSE}
      }
      if(y+1 <= 100 & change){
        if(graph[x,y+1] == 1){change<-FALSE}
      }        
      if(y >= 2 & change){
        if(graph[x,y-1] == 1){change<-FALSE}
      }
      if(change){graph[x,y]<-1}
      nu.i <- sum(graph)^2
      nu.j <- nu.i + runif(1)
      fx[i] <- c(sum(graph))
    }
    else graph[x,y]<-0
    nu.i <- sum(graph)^2
    nu.j <- nu.i + runif(1)
    fx[i] <- c(sum(graph))
  }
  else i = i+1
}
return(fx)
}

x <FUN.COIN1(40000)
plot(x)


# QUESTION 3c
FUN.SAMPLE <- function(step){
  
  #Initailize
  fx <- NULL
  flip <- NULL
  valid <- list()
  for (i in 1:step){
    
    change <- TRUE
    
    #Generate All 10 Samples
    x1 <- sample(x=c(-1,1), size=20, prob=c(.45,.55), replace=TRUE)
    x2 <- sample(x=c(-1,1), size=20, prob=c(.45,.55), replace=TRUE)
    x3 <- sample(x=c(-1,1), size=20, prob=c(.45,.55), replace=TRUE)
    x4 <- sample(x=c(-1,1), size=20, prob=c(.45,.55), replace=TRUE)
    x5 <- sample(x=c(-1,1), size=20, prob=c(.45,.55), replace=TRUE)
    x6 <- sample(x=c(-1,1), size=20, prob=c(.45,.55), replace=TRUE)
    x7 <- sample(x=c(-1,1), size=20, prob=c(.45,.55), replace=TRUE)
    x8 <- sample(x=c(-1,1), size=20, prob=c(.45,.55), replace=TRUE)
    x9 <- sample(x=c(-1,1), size=20, prob=c(.45,.55), replace=TRUE)
    x10 <- sample(x=c(-1,1), size=20, prob=c(.45,.55), replace=TRUE)
    
      # Check 1-20
      if(sum(x1)==2 & change){flip <- x1}
        else{change<-FALSE}
      
      # Check 21-40
      if(sum(x2)==2 & change){flip <- c(x1,x2)}
        else{change<-FALSE}
    
      # Check 41-60
      if(sum(x3)==2 & change){flip <- c(x1,x2,x3)}
        else{change<-FALSE}
    
      # Check 61-80
      if(sum(x4)==2 & change){flip <- c(x1,x2,x3,x4)}
        else{change<-FALSE}
      
      # Check 81-100
      if(sum(x5)==2 & change){flip <- c(x1,x2,x3,x4,x5)}
        else{change<-FALSE}
      
      # Check 101-120
      if(sum(x6)==2 & change){flip <- c(x1,x2,x3,x4,x5,x6)}
        else{change<-FALSE}
      
      # Check 121-140
      if(sum(x7)==2 & change){flip <- c(x1,x2,x3,x4,x5,x6,x7)}
        else{change<-FALSE}
      
      # Check 141-160
      if(sum(x8)==2 & change){flip <- c(x1,x2,x3,x4,x5,x6,x7,x8)}
        else{change<-FALSE}
      
      # Check 161-180
      if(sum(x9)==2 & change){flip <- c(x1,x2,x3,x4,x5,x6,x7,x8,x9)}
        else{change<-FALSE}
      
      # Check 181-200
      if(sum(x10)==2 & change){flip <- c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)}
        else{change<-FALSE}
  
  if(length(valid)==0){append(valid,flip)}
    else{    
      if(flip %in% valid){flip <- NULL}
        else{valid <- append(valid,flip)}
    }
  x <- valid
  }
  return(x)
}