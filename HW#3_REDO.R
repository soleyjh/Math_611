# JAMES HYER SOLEY
# HW #3 STEEPEST ASCENT TEST FUNCTION RE-D0
# 11/09/2014
# ORIGINAL HW DUE IN SEPTEMBER

# Establish the Gradients for Test Function
gradient <- function(x) {
  df.dx1 <- (20 - 2*x[1])
  df.dx2 <- (10 - 2*x[2])
  grad <- c(df.dx1,df.dx2)
  return(grad)
}

# Function to Calculate the gradient Norm
norm <- function(grad) {
  gradnorm = sqrt(sum(grad^2))
  return(gradnorm)
}

##steepest ascent algorithm for test function

testascent <- function(starts = c(0,0), target = .01, step = .1) {
  x.new <- starts
  grad.new <- gradient(x.new)
  path <- matrix(0,1,2)
  gradients <- matrix(0,1,2)
  
  while (norm(grad.new) > target) {
    
    if (grad.new[1] %in% gradients[,1] && grad.new[2] %in% gradients[,2]) {
      step <- step/2
    }
    
    path <- rbind(path, x.new)  #add new x to previous path
    
    x.old <- x.new 
    grad.old <- grad.new 
    gradients <- rbind(gradients,grad.old)
    
    dir <- grad.old/norm(grad.old)
    
    x.new <- x.old + dir*step
    grad.new <- gradient(x.new)
  }
  
  return(rbind(path,x.new))
}

test <- testascent()