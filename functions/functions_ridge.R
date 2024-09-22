#### Ridge Regression

#standardize data

standardize_data <- function(X){
  
  return(as.matrix(as.data.frame(scale(X))))
}

ridge_regression <- function(Y, X, mu, N = 272){
  
  
  beta <- solve(t(X)%*%X + mu*diag(N)) %*% t(X) %*% Y
  
  return(beta)
  
}


