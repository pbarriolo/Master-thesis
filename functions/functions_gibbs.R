# Functions for gibbs sampler of Korobolis Paper
###############

# initializing functions

initiate_gamma <- function(beta_0) {
  gamma_ <- rep(0, length(beta_0))

  for (i in 1:length(gamma_)) {
    if (beta_0[i] <= 0.0001) {
      gamma_[i] <- 0
    } else {
      gamma_[i] <- 1
    }
  }
  return(gamma_)
}


initialize_gibbs <- function(Y, X, T, N, a_1, a_2, b_1, b_2) {
  mu <- 0.5

  beta_0 <- rnorm(ncol(X), mean = 0, sd = .25) # rep(0, ncol(X))  all betas equal to zero

  var_residuals <- t(Y[1:T] - X[1:T, ] %*% beta_0) %*% (Y[1:T] - X[1:T, ] %*% beta_0) / (T - N - 1)

  delta_0 <- (1 / rgamma(N, a_1, rate = a_2)) # diag(var_residuals[1] * solve(t(X[1:T, ]) %*% X[1:T, ] + mu * diag(N)) %*% t(X[1:T, ]) %*% X[1:T, ] %*% solve(t(X[1:T, ]) %*% X[1:T, ] + mu * diag(N)))

  gamma_0 <- initiate_gamma(beta_0)

  Z_0 <- rep(1, T) # rexp(T)*0

  pi_0 <- rbeta(1, b_1, b_2)

  return(list(beta_0, pi_0, delta_0, gamma_0, Z_0))
}

# drawing functions

draw_beta <- function(Y, X, Z, thau_2, theta, delta, gamma_, c_, T) {
  new_delta <- gamma_ * delta + c_ * (1 - gamma_) * delta

  # V_inv <- solve(diag(new_delta))

  # temp_mean <- 0

  # for (t in 1:T) {
  #  V_inv <- V_inv + X[t, ] %*% t(X[t, ]) / (thau_2 * Z[t])
  #  temp_mean <- temp_mean + X[t, ] * (Y[t] - theta * Z[t]) / (thau_2 * Z[t])
  # }
  V_inv <- solve(diag(new_delta), tol = 1e-20) + t(X) %*% diag(1 / (thau_2 * Z)) %*% X
  temp_mean <- t(X) %*% diag(1 / (thau_2 * Z)) %*% (Y - theta * Z)

  V <- solve(V_inv + diag(rep(1e-20, ncol(X))), tol = 1e-20)
  mean <- V %*% temp_mean

  # print(V)
  # print(mean)

  return(rmvnorm(1, mean, V))
}

draw_delta_i <- function(beta_i, a_1, a_2) {
  new_a_1 <- a_1 + 0.5
  new_a_2 <- a_2 + 0.5 * beta_i^2

  return(1 / rgamma(1, new_a_1, rate = new_a_2))
}

draw_delta <- function(beta, a_1, a_2) {
  sapply(beta, function(x) draw_delta_i(x, a_2, a_2))
}

compute_likelihood <- function(Y, X, Z, beta, thau_2, theta, T, n) {
  # temp_1 <- 1
  # temp_2 <- 0
  # for (t in 1:T) {
  #  temp_1 <- temp_1 * (Z[t])^(-0.5)
  #  temp_2 <- temp_2 + (Y[t] - theta * Z[t] - t(X[t, ]) %*% t(beta))^2 / (thau_2 * Z[t])
  # }

  temp_1 <- prod(Z^(-0.5))

  temp_2 <- sum((Y - theta * Z - X %*% t(beta))^2 / (thau_2 * Z))

  # cat("
  # temp1 of likelihood", temp_1)
  # cat("
  # temp2 of likelihood", temp_2)

  return(temp_1 * exp(-0.5 * temp_2))
}

draw_gamma_ <- function(Y, X, Z, pi_0, thau_2, theta, delta, gamma_, c_, T, n) {
  gamma_[1] <- 0
  beta <- draw_beta(Y, X, Z, thau_2, theta, delta, gamma_, c_, T)
  likelihood_0 <- compute_likelihood(Y, X, Z, beta, thau_2, theta, T, n)

  gamma_[1] <- 1
  beta <- draw_beta(Y, X, Z, thau_2, theta, delta, gamma_, c_, T)
  likelihood_1 <- compute_likelihood(Y, X, Z, beta, thau_2, theta, T, n)

  temp_1 <- log(pi_0 * likelihood_1)
  temp_2 <- log(pi_0 * likelihood_1 + (1 - pi_0) * likelihood_0)

  # cat("
  # element of gamma_: 1")
  # cat("
  # Likelihood_1", likelihood_0)
  # cat("
  # Likelihood_0", likelihood_1)


  # cat("
  #    Intial temps",temp_1 , temp_2)

  pi_tilde <- exp(temp_1 - temp_2)

  gamma_[1] <- rbernoulli(1, pi_tilde)

  # cat("
  #  gamma_:", gamma_)

  for (i in 2:n) {
    # likelihoods <- c()
    # for (j in 1:2) {
    #  gamma_[i] <- j - 1
    #  beta <- draw_beta(Y, X, Z, thau_2, theta, delta, gamma_, c_, T)

    # likelihood <- compute_likelihood(Y, X, Z, beta, thau_2, theta, T, n)
    #  likelihoods <- append(likelihoods, likelihood)
    # }
    if (gamma_[i] == 0) {
      if (gamma_[i - 1] == 0) {
        likelihood_0 <- likelihood_0
      } else {
        likelihood_0 <- likelihood_1
      }

      gamma_[i] <- 1
      beta <- draw_beta(Y, X, Z, thau_2, theta, delta, gamma_, c_, T)
      likelihood_1 <- compute_likelihood(Y, X, Z, beta, thau_2, theta, T, n)
    } else {
      if (gamma_[i - 1] == 0) {
        likelihood_1 <- likelihood_0
      } else {
        likelihood_1 <- likelihood_1
      }

      gamma_[i] <- 0
      beta <- draw_beta(Y, X, Z, thau_2, theta, delta, gamma_, c_, T)
      likelihood_0 <- compute_likelihood(Y, X, Z, beta, thau_2, theta, T, n)
    }


    temp_1 <- log(pi_0 * likelihood_1)
    temp_2 <- log(pi_0 * likelihood_1 + (1 - pi_0) * likelihood_0)

    # cat("
    # element of gamma_", i)
    # cat("
    # Likelihood_1", likelihood_0)
    # cat("
    # Likelihood_0", likelihood_1)

    pi_tilde <- exp(temp_1 - temp_2)

    gamma_[i] <- rbernoulli(1, pi_tilde)
    # print(i)
    # print(gamma_)
  }
  return(list(gamma_, likelihood_0, likelihood_1))
}

draw_gamma_i_new <- function(Y, X, Z, pi_0, thau_2, theta, delta, gamma_, c_, T, n, i) {
  gamma_[i] <- 0
  beta <- draw_beta(Y, X, Z, thau_2, theta, delta, gamma_, c_, T)
  likelihood_0 <- compute_likelihood(Y, X, Z, beta, thau_2, theta, T, n)

  gamma_[i] <- 1
  beta <- draw_beta(Y, X, Z, thau_2, theta, delta, gamma_, c_, T)
  likelihood_1 <- compute_likelihood(Y, X, Z, beta, thau_2, theta, T, n)

  temp_1 <- pi_0 * likelihood_1
  temp_2 <- pi_0 * likelihood_1 + (1 - pi_0) * likelihood_0

  pi_tilde <- temp_1 / temp_2

  gamma_[i] <- rbernoulli(1, pi_tilde)

  return(gamma_)
}

# Function to sample gamma_p (corrected)
sample_gamma_new_v2 <- function(beta_p, Y, X, Z, pi_0, thau_2, theta, delta, gamma_p, c_, T, n, current_index) {
  # Set gamma_p[current_index] to 1 (included) and calculate likelihood
  gamma_p_1 <- gamma_p
  gamma_p_1[current_index] <- 1
  beta_p_1 <- beta_p
  beta_p_1[current_index] <- beta_p[current_index] # Include beta_p[current_index]

  likelihood_gamma1 <- exp(-sum((Y - X %*% t(beta_p_1) - theta * Z)^2 / (2 * thau_2 * Z)))

  # Set gamma_p[current_index] to 0 (excluded) and calculate likelihood
  gamma_p_0 <- gamma_p
  gamma_p_0[current_index] <- 0
  beta_p_0 <- beta_p
  beta_p_0[current_index] <- beta_p[current_index] * c_ # Exclude beta_p[current_index]

  likelihood_gamma0 <- exp(-sum((Y - X %*% t(beta_p_0) - theta * Z)^2 / (2 * thau_2 * Z)))

  # Calculate the posterior inclusion probability for gamma_p[current_index]
  pi_gamma_p <- pi_0 * likelihood_gamma1 / (pi_0 * likelihood_gamma1 + (1 - pi_0) * likelihood_gamma0)

  # Sample gamma_p[current_index] from Bernoulli(pi_gamma_p)
  return(list(rbinom(1, 1, pi_gamma_p), likelihood_gamma0, likelihood_gamma1))
}



draw_pi <- function(gamma_, b_1, b_2, n) {
  n_gamma_ <- gamma_ %*% rep(1, n)

  new_b_1 <- b_1 + n_gamma_
  new_b_2 <- b_2 + n - n_gamma_

  return(rbeta(1, new_b_1, new_b_2))
}

draw_Z_og <- function(Y, X, beta, thau_2, theta, T) {
  # NO FUNCIONA!!!
  # for (t in 1:T) {
  #  kappa_1 <- kappa_1 + (Y[t] - t(X[t, ]) %*% t(beta)) / sqrt(thau_2)
  # }
  kappa_1 <- (colSums(Y - X %*% t(beta)) / sqrt(thau_2))

  kappa_2 <- sqrt((2 + theta^2) / thau_2)

  return(rgig(T, 0.5, kappa_1^2, kappa_2^2))
}

draw_new_Z <- function(Y, X, beta, thau_2, theta, T) {
  Z_temp <- rep(1, T)

  kappa_2 <- sqrt((2 + theta^2) / thau_2) # sqrt((2*thau_2 + theta^2) / thau_2)

  for (t in 1:T) {
    kappa_1 <- abs(Y[t] - X[t, ] %*% t(beta)) / sqrt(thau_2)

    # cat("kappas:", kappa_1,",", kappa_2)

    Z_temp[t] <- rgig(1, 0.5, kappa_1^2, kappa_2^2)
    # cat("
    # Z at t:", t, " is ", Z_temp[t])
  }


  return(Z_temp)
}






# gibbs sampler function

gibbs_sampling <- function(Y, X, Z, pi_0, thau_2, theta, beta, gamma_, delta, T, n, iterations, burn, a_1, a_2, b_1, b_2, c_) {
  beta_list <- data.frame(matrix(nrow = iterations + burn, ncol = ncol(X)))
  gamma_list <- data.frame(matrix(nrow = iterations + burn, ncol = ncol(X)))
  delta_list <- data.frame(matrix(nrow = iterations + burn, ncol = ncol(X)))
  Z_list <- data.frame(matrix(nrow = iterations + burn, ncol = T))
  pi_0_list <- data.frame(matrix(nrow = iterations + burn, ncol = 1))
  l_0_list <- data.frame(matrix(nrow = iterations + burn, ncol = 1))
  l_1_list <- data.frame(matrix(nrow = iterations + burn, ncol = 1))


  for (i in 1:(iterations + burn)) {
    if (i %% 100 == 0) {
      cat("
        -----------------------------ITERATION:", i, "
        ")
    }

    beta <- draw_beta(Y[1:T], X[1:T, ], Z[1:T], thau_2, theta, delta, gamma_, c_, T)


    delta <- draw_delta(beta, a_1, a_2)

    # results_gamma_ <- draw_gamma_(Y, X[1:T, ], Z[1:T], pi_0, thau_2, theta, delta, gamma_, c_, T, n)

    # gamma_ <- results_gamma_[[1]]

    for (j in 1:ncol(X)) {
      results_gamma_ <- sample_gamma_new_v2(beta, Y[1:T], X[1:T, ], Z[1:T], pi_0, thau_2, theta, delta, gamma_, c_, T, n, j)
      gamma_[j] <- results_gamma_[[1]]
    }

    pi_0 <- draw_pi(gamma_, b_1, b_2, n)

    Z <- draw_new_Z(Y[1:T], X[1:T, ], beta, thau_2, theta, T)

    # if (i > burn) {
    beta_list[i, ] <- beta
    gamma_list[i, ] <- gamma_
    delta_list[i, ] <- delta
    Z_list[i, ] <- Z
    pi_0_list[i, ] <- pi_0
    l_0_list[i, ] <- results_gamma_[[2]]
    l_1_list[i, ] <- results_gamma_[[3]]

    # if ((i %% 10000 == 0) | i == 2 | i == 100 | i == 500 | i == 1000 | i == 5000) {
    #  write.csv(beta_list, paste("results/tests/", vintage_date, "_h", horizon, "_beta_p_", p_level, "_u_u_og.csv", sep = ""))
    #  write.csv(gamma_list, paste("results/tests/", vintage_date, "_h", horizon, "_gamma_p_", p_level, "_u_u_og.csv", sep = ""))
    #  write.csv(delta_list, paste("results/tests/", vintage_date, "_h", horizon, "_delta_p_", p_level, "_u_u_og.csv", sep = ""))
    #  write.csv(Z_list, paste("results/tests/", vintage_date, "_h", horizon, "_Z_p_", p_level, "_u_u_og.csv", sep = ""))
    #  write.csv(pi_0_list, paste("results/tests/", vintage_date, "_h", horizon, "_pi0_p_", p_level, "_u_u_og.csv", sep = ""))
    #  write.csv(l_0_list, paste("results/tests/", vintage_date, "_h", horizon, "_l0_p_", p_level, "_u_u_og.csv", sep = ""))
    #  write.csv(l_1_list, paste("results/tests/", vintage_date, "_h", horizon, "_l1_p_", p_level, "_u_u_og.csv", sep = ""))
    # }
  }

  return(list(beta_list, gamma_list, delta_list, Z_list, pi_0_list, l_0_list, l_1_list))
}



# gibbs sampler function for prediction

gibbs_sampling_prediction <- function(Y, X, Z, pi_0, thau_2, theta, beta, gamma_, delta, T, n, iterations, burn, a_1, a_2, b_1, b_2, c_) {
  # we make an empty matrix to store results
  beta_list <- data.frame(matrix(nrow = iterations, ncol = ncol(X)))
  gamma_list <- data.frame(matrix(nrow = iterations, ncol = ncol(X)))
  delta_list <- data.frame(matrix(nrow = iterations, ncol = ncol(X)))

  for (i in 1:(iterations + burn)) {
    if (i %% 100 == 0) {
      cat("
        -----------------------------ITERATION:", i, "
        ")
    }

    delta <- draw_delta(beta, a_1, a_2)

    gamma_ <- draw_gamma_(Y[1:T], X[1:T, ], Z[1:T], pi_0, thau_2, theta, delta, gamma_, c_, T, n)

    beta <- draw_beta(Y[1:T], X[1:T, ], Z[1:T], thau_2, theta, delta, gamma_, c_, T)

    pi_0 <- draw_pi(gamma_, b_1, b_2, n)

    Z <- draw_new_Z(Y[1:T], X[1:T, ], beta, thau_2, theta, T)

    beta_list[i, ] <- beta
    gamma_list[i, ] <- gamma_
    delta_list[i, ] <- delta
  }

  return(list(beta_list, gamma_list, delta_list))
}

gibbs_sampling_old <- function(Y, X, Z, pi_0, thau_2, theta, beta, gamma_, delta, T, n, iterations, burn, a_1, a_2, b_1, b_2, c_) {
  # we make an empty matrix to store results
  beta_list <- data.frame(matrix(nrow = iterations, ncol = ncol(X)))
  gamma_list <- data.frame(matrix(nrow = iterations, ncol = ncol(X)))
  delta_list <- data.frame(matrix(nrow = iterations, ncol = ncol(X)))

  for (i in 1:(iterations + burn)) {
    if (i %% 100 == 0) {
      cat("
        -----------------------------  ITERATION:", i, " -------------
        ")
    }

    beta <- draw_beta(Y[1:T], X[1:T, ], Z[1:T], thau_2, theta, delta, gamma_, c_, T)

    delta <- draw_delta(beta, a_1, a_2)

    gamma_ <- draw_gamma_(Y[1:T], X[1:T, ], Z[1:T], pi_0, thau_2, theta, delta, gamma_, c_, T, n)

    pi_0 <- draw_pi(gamma_, b_1, b_2, n)

    Z <- draw_new_Z(Y[1:T], X[1:T, ], beta, thau_2, theta, T)

    beta_list[i, ] <- beta
    gamma_list[i, ] <- gamma_
    delta_list[i, ] <- delta
  }

  return(list(beta_list, gamma_list, delta_list))
}

gibbs_sampling_prediction_v2 <- function(Y, X, Z, pi_0, thau_2, theta, beta, gamma_, delta, T, n, iterations, burn, a_1, a_2, b_1, b_2, c_) {
  # we make an empty matrix to store results
  beta_list <- data.frame(matrix(nrow = iterations, ncol = ncol(X)))
  gamma_list <- data.frame(matrix(nrow = iterations, ncol = ncol(X)))
  delta_list <- data.frame(matrix(nrow = iterations, ncol = ncol(X)))

  for (i in seq(1, (iterations + burn), ncol(X))) {
    if (i %% 100 == 0) {
      cat("
        ------------------------------- ITERATION:", i, " -------------
        ")
    }
    for (j in 1:ncol(X)) {
      beta <- draw_beta(Y[1:T], X[1:T, ], Z[1:T], thau_2, theta, delta, gamma_, c_, T)

      delta <- draw_delta(beta, a_1, a_2)

      gamma_ <- draw_gamma_i_new(Y[1:T], X[1:T, ], Z[1:T], pi_0, thau_2, theta, delta, gamma_, c_, T, n, j)

      pi_0 <- draw_pi(gamma_, b_1, b_2, n)

      Z <- draw_new_Z(Y[1:T], X[1:T, ], beta, thau_2, theta, T)

      beta_list[i + j - 1, ] <- beta
      gamma_list[i + j - 1, ] <- gamma_
      delta_list[i + j - 1, ] <- delta
    }
  }



  return(list(beta_list, gamma_list, delta_list))
}