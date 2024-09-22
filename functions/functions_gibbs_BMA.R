# Functions for gibbs sampler for BMA for a simple mean regression
###############

# initializing functions

initiate_gamma_bma <- function(beta_0) {
    gamma_ <- rep(0, length(beta_0))

    for (i in 1:length(gamma_)) {
        if (beta_0[i] <= 0.001) {
            gamma_[i] <- 0
        } else {
            gamma_[i] <- 1
        }
    }
    return(gamma_)
}


initialize_gibbs_bma <- function(Y, X, T, N, b_1, b_2) {
    mu <- 0.5

    beta_0 <- rep(0, ncol(X)) # all betas equal to zero

    var_residuals <- t(Y[1:T] - X[1:T, ] %*% beta_0) %*% (Y[1:T] - X[1:T, ] %*% beta_0) / (T - N - 1)

    delta_0 <- diag(var_residuals[1] * solve(t(X[1:T, ]) %*% X[1:T, ] + mu * diag(N)) %*% t(X[1:T, ]) %*% X[1:T, ] %*% solve(t(X[1:T, ]) %*% X[1:T, ] + mu * diag(N)))

    gamma_0 <- initiate_gamma(beta_0)

    pi_0 <- 0.5

    return(list(beta_0, delta_0, gamma_0, pi_0))
}

# draw_bmaing functions

draw_beta_bma <- function(Y, X, sigma_2, delta, gamma_, c_, T) {
    # cat("
    # delta:", delta)
    new_delta <- gamma_ * delta + c_ * (1 - gamma_) * delta
    # cat("
    # new_delta:", new_delta)

    # V_inv <- solve(diag(new_delta))

    # temp_mean <- 0

    # for (t in 1:T) {
    #  V_inv <- V_inv + X[t, ] %*% t(X[t, ]) / (thau_2 * Z[t])
    #  temp_mean <- temp_mean + X[t, ] * (Y[t] - theta * Z[t]) / (thau_2 * Z[t])
    # }
    V_inv <- solve(diag(new_delta), tol = 1e-20) + t(X) %*% diag(1 / (sigma_2), T) %*% X
    temp_mean <- t(X) %*% ((Y) / (sigma_2))

    V <- solve(V_inv, tol = 1e-20)
    mean <- V %*% temp_mean

    # print(V)
    # print(mean)

    return(rmvnorm(1, mean, V))
}

draw_delta_i_bma <- function(beta_i, a_1, a_2) {
    new_a_1 <- a_1 + 0.5
    new_a_2 <- a_2 + 0.5 * beta_i^2

    return(1 / rgamma(1, new_a_1, rate = new_a_2))
}

draw_delta_bma <- function(beta, a_1, a_2) {
    sapply(beta, function(x) draw_delta_i(x, a_2, a_2))
}

compute_likelihood_bma <- function(Y, X, beta, sigma_2, T, n) {
    # temp_1 <- 1
    # temp_2 <- 0
    # for (t in 1:T) {
    #  temp_1 <- temp_1 * (Z[t])^(-0.5)
    #  temp_2 <- temp_2 + (Y[t] - theta * Z[t] - t(X[t, ]) %*% t(beta))^2 / (thau_2 * Z[t])
    # }

    temp_2 <- sum((Y - X %*% t(beta))^2 / (sigma_2))

    return(exp(-0.5 * temp_2))
}

draw_gamma_bma <- function(Y, X, pi_0, sigma_2, delta, gamma_, c_, T, n) {
    gamma_[1] <- 0
    beta <- draw_beta_bma(Y, X, sigma_2, delta, gamma_, c_, T)
    likelihood_0 <- compute_likelihood_bma(Y, X, beta, sigma_2, T, n)

    gamma_[1] <- 1
    beta <- draw_beta_bma(Y, X, sigma_2, delta, gamma_, c_, T)
    likelihood_1 <- compute_likelihood_bma(Y, X, beta, sigma_2, T, n)

    temp_1 <- pi_0 * likelihood_1
    temp_2 <- pi_0 * likelihood_1 + (1 - pi_0) * likelihood_0

    # cat("
    #    Intial temps",temp_1 , temp_2)

    pi_tilde <- temp_1 / temp_2

    gamma_[1] <- rbernoulli(1, pi_tilde)

    for (i in 2:n) {
        # likelihoods <- c()
        # for (j in 1:2) {
        #  gamma_[i] <- j - 1
        #  beta <- draw_beta(Y, X, sigma_2, delta, gamma_, c_, T)

        # likelihood <- compute_likelihood(Y, X, beta, sigma_2, T, n)
        #  likelihoods <- append(likelihoods, likelihood)
        # }
        if (gamma_[i] == 0) {
            if (gamma_[i - 1] == 0) {
                likelihood_0 <- likelihood_0
            } else {
                likelihood_0 <- likelihood_1
            }

            gamma_[i] <- 1
            beta <- draw_beta_bma(Y, X, sigma_2, delta, gamma_, c_, T)
            likelihood_1 <- compute_likelihood_bma(Y, X, beta, sigma_2, T, n)
        } else {
            if (gamma_[i - 1] == 0) {
                likelihood_1 <- likelihood_0
            } else {
                likelihood_1 <- likelihood_1
            }

            gamma_[i] <- 0
            beta <- draw_beta_bma(Y, X, sigma_2, delta, gamma_, c_, T)
            likelihood_0 <- compute_likelihood_bma(Y, X, beta, sigma_2, T, n)
        }


        temp_1 <- pi_0 * likelihood_1
        temp_2 <- pi_0 * likelihood_1 + (1 - pi_0) * likelihood_0
        # cat("element of gamma_", i)
        # cat(" Likelihood_1", likelihood_0)
        # cat(" Likelihood_0", likelihood_1)
        pi_tilde <- temp_1 / temp_2

        gamma_[i] <- rbernoulli(1, pi_tilde)
        # print(i)
        # print(gamma_)
    }
    return(gamma_)
}

draw_pi_bma <- function(gamma_, b_1, b_2, n) {
    n_gamma_ <- gamma_ %*% rep(1, n)

    new_b_1 <- b_1 + n_gamma_
    new_b_2 <- b_2 + n - n_gamma_

    return(rbeta(1, new_b_1, new_b_2))
}

draw_forecast_Y_new_bma <- function(X_out, beta, sigma_2) {
    u <- rnorm(1, mean = 0, sd = 1)

    return(X_out %*% t(beta) + sqrt(sigma_2) * u)
}

sample_gamma_new_v2_bma <- function(beta_p, Y, X, pi_0, delta, gamma_p, c_, T, n, current_index) {
    # Set gamma_p[current_index] to 1 (included) and calculate likelihood
    gamma_p_1 <- gamma_p
    gamma_p_1[current_index] <- 1
    beta_p_1 <- beta_p
    beta_p_1[current_index] <- beta_p[current_index] # Include beta_p[current_index]


    likelihood_gamma1 <- exp(-sum((Y - X %*% t(beta_p_1))^2 / 2))

    # Set gamma_p[current_index] to 0 (excluded) and calculate likelihood
    gamma_p_0 <- gamma_p
    gamma_p_0[current_index] <- 0
    beta_p_0 <- beta_p
    beta_p_0[current_index] <- beta_p[current_index] * c_ # Exclude beta_p[current_index]

    likelihood_gamma0 <- exp(-sum((Y - X %*% t(beta_p_0))^2 / 2))
    cat("
    li0:", likelihood_gamma0)
    cat("
    li1:", likelihood_gamma1)

    # Calculate the posterior inclusion probability for gamma_p[current_index]
    pi_gamma_p <- pi_0 * likelihood_gamma1 / (pi_0 * likelihood_gamma1 + (1 - pi_0) * likelihood_gamma0)

    # Sample gamma_p[current_index] from Bernoulli(pi_gamma_p)
    return(rbinom(1, 1, pi_gamma_p))
}

# gibbs samp_bmaler function for prediction

gibbs_sampling_bma <- function(Y, X, sigma_2, beta, delta, gamma_, pi_0, T, n, iterations, burn, a_1, a_2, b_1, b_2, c_) {
    # we make an empty matrix to store results
    beta_list <- data.frame(matrix(nrow = iterations, ncol = ncol(X)))
    gamma_list <- data.frame(matrix(nrow = iterations, ncol = ncol(X)))
    Y_list <- data.frame(matrix(nrow = iterations, ncol = 1))

    for (i in 1:(iterations + burn)) {
        if (i %% 100 == 0) {
            cat("
        -----------------------------ITERATION:", i, "
        ")
        }

        beta <- draw_beta_bma(Y[1:T], X[1:T, ], sigma_2, delta, gamma_, c_, T)

        delta <- draw_delta_bma(beta, a_1, a_2)

        gamma_ <- draw_gamma_bma(Y[1:T], X[1:T, ], pi_0, sigma_2, delta, gamma_, c_, T, n)

        #    for (j in 1:ncol(X)) {
        #        gamma_[j] <- sample_gamma_new_v2_bma(beta, Y, X[1:T, ], pi_0, delta, gamma_, c_, T, n, j)
        #    }
        #    cat("
        # gamma:", gamma_)

        # pi_0 <- draw_pi_bma(gamma_, b_1, b_2, n)

        Y_h <- draw_forecast_Y_new_bma(X[T + h + 1, ], beta, sigma_2)

        beta_list[i - burn, ] <- beta
        gamma_list[i - burn, ] <- gamma_
        Y_list[i, ] <- Y_h
    }

    return(list(beta_list, gamma_list, Y_list))
}