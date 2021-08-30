library(tidyverse)

Newton_Raphson_Poisson <- function(beta_0, data, epsilon){
  # set initial beta guess
  beta <- beta_0
  
  # get the predictor values
  X <- data %>% 
    select(2:5) %>% 
    as.matrix()
  
  # get the response variable (crashes) values
  y <- data %>% 
    select(6) %>% 
    as.matrix()
  
  # convergence flag
  flag <- 0
  
  while(flag != 1){
    # save the calc.score value
    score <- calc.score(beta, X, y)
    
    # get the new beta value according to Newton Raphson method
    beta_prop <- beta - solve(calc.hess(beta, X)) %*% score
    
    #check for convergence
    if(all(score[1] < epsilon, score[2] < epsilon, score[3] < epsilon, 
           score[4] < epsilon)){
      flag <- 1
    }
    else{
      # update the new beta value
      beta <- beta_prop
    }
  }
  
  return(beta)
}

calc.hess <- function(beta, X){
  d1 <- rep(0, 9) %>% 
    matrix(nrow = 3)
  
  for(i in 1: dim(X)[[1]]){
    d1 <- d1 - exp(X[i,] %*% beta) %*% X[i,] %*% t(X[i,])
  }
  return(d1)
}

calc.score <- function(beta, X, y){
  d2 <- rep(0, length(beta))
  for(i in 1:length(y)){
    d2 <- d2 + (y[i] - exp(t(X[i,]) %*% beta)) %*% X[i,]
  }
  return(colSums(d2))
}
