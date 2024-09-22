###### functions for simulating data for the QR-BMA

simulate_data <- function(T, k, s, rho, sigma = 1) {
    # Step 1: Set up parameters
    set.seed(123) # For reproducibility

    # Step 2: Generate the Toeplitz correlation matrix
    # This creates a k x k matrix with corr(x_t,i, x_t,j) = rho^|i-j|
    toeplitz_matrix <- toeplitz(rho^(0:(k - 1)))

    # Step 3: Simulate the predictors x_t
    X <- mvrnorm(n = T, mu = rep(0, k), Sigma = toeplitz_matrix)

    # Step 4: Generate the coefficients beta
    beta <- rep(0, k) # Start with a beta vector of zeros
    non_zero_indices <- sample(1:k, s) # Randomly select s indices
    beta[non_zero_indices] <- rnorm(s, mean = 0, sd = 1) # Set these to N(0,1)

    # Step 5: Simulate the error term epsilon
    epsilon <- rnorm(T, mean = 0, sd = sqrt(sigma))

    # Step 6: Simulate the dependent variable y_t
    y <- X %*% beta + epsilon

    return(list(y, X, beta))
}

simulate_data_one_value <- function(T, k, s, rho, sigma, iterations = 1000) {
    simulations <- simulate_data(T, k, s, rho, sigma)

    y <- simulations[[1]]
    X <- simulations[[2]]
    beta <- simulations[[3]]

    # Step 1: Set up parameters
    set.seed(123) # For reproducibility

    element <- rdunif(1, 1, T)
    x_it <- matrix(rep(X[element, ], each = iterations), nrow = iterations, byrow = FALSE)

    # Step 4: Generate the coefficients beta
    beta <- rep(0, k) # Start with a beta vector of zeros
    non_zero_indices <- sample(1:k, s) # Randomly select s indices
    beta[non_zero_indices] <- rnorm(s, mean = 0, sd = 1) # Set these to N(0,1)

    # Step 5: Simulate the error term epsilon
    epsilon_it <- rnorm(iterations, mean = 0, sd = sqrt(sigma))

    # Step 6: Simulate the dependent variable y_t
    y_it <- x_it %*% beta + epsilon_it

    return(list(y, X, beta, y_it, x_it[1, 1:ncol(x_it)]))
}