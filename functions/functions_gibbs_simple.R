###############Functions for gibbs sampler of simplified Korobolis Paper
############### gamma = 1 a.s.
########### having gamma = 1 a.s. means that we ignore two posterior distributions: for gamma and for pi_0
#initiate functions

initiate_gamma <- function(beta_0){
  gamma = rep(0,length(beta_0))
  for ( i in 1:length(gamma)){
    gamma[i] <- 1
  }
  return (gamma)
}

#drawing functions

draw_beta <- function(Y, X, Z, thau_2, delta, gamma, c_, T){
  
  V_inv <- solve(diag(delta))
  
  temp_mean = 0
  
  for (t in 1:T){
    V_inv = V_inv + X[t,]%*%t(X[t,])/(thau_2*Z[t])
    temp_mean = temp_mean + X[t,]*(Y[t] - thau_2*Z[t])/(thau_2*Z[t])
  }
  #print(det(V_inv))
  V = solve(V_inv)
  mean = V%*%temp_mean
  
  return (rmvnorm(1, mean, V))
}

draw_delta_i <- function(beta_i, a_1, a_2){
  
  new_a_1 = a_1 + 0.5
  new_a_2 = a_2 + 0.5*beta_i^2
  
  return( 1/rgamma(1, new_a_1, new_a_2))
  
}

draw_delta <- function(beta, a_1, a_2){
  
  sapply(beta, function(x) draw_delta_i(x, a_2, a_2))
  
}


draw_Z <- function(Y, X, beta, thau_2, theta, T){
  
  kappa_1 = 0
  
  for (t in 1:T){
    kappa_1 = kappa_1 + abs(Y[t] - t(X[t,])%*%t(beta))/sqrt(thau_2)
  }
  
  kappa_2 = sqrt((2+theta^2)/thau_2)
  
  cat("kappas:", abs(kappa_1),",", kappa_2)
  
  return (rgig(T, 0.5, kappa_1, kappa_2))
  
}

draw_new_Z <- function(Y, X, beta, thau_2, theta, T){
  
  Z_temp = rep(1,T)
  
  kappa_2 = sqrt((2+theta^2)/thau_2)
  
  for (t in 1:T){
    
    kappa_1 =abs(Y[t] - t(X[t,])%*%t(beta))/sqrt(thau_2)
    
    #cat("kappas:", kappa_1,",", kappa_2)
    
    Z_temp[t] = rgig(1, 0.5, kappa_1, kappa_2)
    
  }
  
  return (Z_temp)
  
}

draw_forecast <- function(X, Z, beta, thau_2, theta, T){
  
  return (rnorm(1, mean = t(X[T,])%*%t(beta) + theta*Z[T], sd = thau_2*Z[T]))
}




# gibbs sampler function

gibbs_sampling_simplified <- function(Y, X, Z, pi_0, thau_2, theta, gamma, delta, T, n, iterations, burn, a_1, a_2, b_1, b_2, c_){
  
  beta_list= data.frame(matrix(nrow = iterations, ncol = ncol(X)))
  Y_h_list= data.frame(matrix(nrow = iterations, ncol = 1))
  
  for(i in 1:iterations + burn){
    cat("
        --------------------ITERATION n:", i)
    
    beta = draw_beta(Y, X, Z, thau_2, delta, gamma, c_, T)
    
    delta = draw_delta(beta, a_1, a_2)
    
    #cat(" delta", delta)
    
    #gamma = draw_gamma(Y, X, Z, pi_0, thau_2, theta, delta, gamma, c_, T, n)
    
    #cat(" gamma:", gamma)
    
    #pi_0 = draw_pi(gamma, b_1, b_2, n)
    
    #cat(" pi_0", pi_0)
    
    Z = draw_new_Z(Y, X, beta, thau_2, theta, T)
    
    #cat(" Z:", Z)
    
    if(i > burn){
      
      beta_list[i-burn,] <-- beta
      Y_h_list[i-burn,] <-- draw_forecast(X, Z, beta, thau_2, theta, T)
      
    }
  }
  
  colnames(beta_list) = colnames(X)
  
  return (list(beta_list, Y_h_list))
}
