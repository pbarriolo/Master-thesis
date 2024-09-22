##########################################################################
################### Run QR-BMA on Simulated data #############################3


# In this file we run the QR-BMA on simulated data
###################################################
library(parallel)
library(foreach)
library(doParallel)
library(readxl)
library(zoo)
library(tidyverse)
library(mvtnorm)
library(ghyp)
library(robustbase)
library(dplyr)

# sources for function
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_transform.R")
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_gibbs.R")
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_import_data.R")
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_ridge.R")

# Import df into one df_frame
path <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/"
setwd(path)
TT <- 50

Y <- as.matrix(read.csv("data/simulations/simulations_T_50_s_5_rho_0.65_sigma_5/y.csv")[, 2])
X <- as.matrix(read.csv("data/simulations/simulations_T_50_s_5_rho_0.65_sigma_5/X.csv")[, 2:31])
x_it <- as.matrix(read.csv("data/simulations/simulations_T_50_s_5_rho_0.65_sigma_5/x_it.csv")[, 2])
Y_it <- as.matrix(read.csv("data/simulations/simulations_T_50_s_5_rho_0.65_sigma_5/y_it.csv")[, 2])
beta_true <- read.csv("data/simulations/simulations_T_50_s_5_rho_0.65_sigma_5/beta_true.csv")[, 2]
N <- ncol(X)

numCores <- detectCores()
registerDoParallel(3)

############################## Initialization  ###################################
# the benchmarck is uninformative priors

# Variables and Parameters definition
ITER <- 2000
BURN <- 500

C_ <- 0.00001 # try with 0
A_1 <- 1 # test with high variance (1, 10) and medium variance (2, 1), (2, 0.1) low variance, uninformative (0.1,0.1)
A_2 <- 10
B_1 <- 1 # ?? test with (5,1) loose, (1,5) restrictive, and (1, 20) very restrictive, uninformative (1,1)
B_2 <- 20 # ??

pi_0 <- rbeta(1, B_1, B_2)

T_list <- (45:49) # (as.integer(length(Y) - length(Y) / 10):as.integer(length(Y) - 1))
p_list <- (1:9) * 0.1

# param_list <- list(c(46, 0.2), c(46, 0.3), c(47, 0.8))

param_list <- expand.grid(T = T_list, p = p_list)

foreach(i = 1:nrow(param_list), .combine = c, .packages = c("readxl", "zoo", "tidyverse", "mvtnorm", "ghyp", "robustbase", "dplyr")) %dopar% {
    T <- param_list[i, "T"]
    p <- param_list[i, "p"]
    # T <- param_list[[i]][1]
    # p <- param_list[[i]][2]

    initial_values <- initialize_gibbs(Y, X, T, N, A_1, A_2, B_1, B_2)

    THETA <- (1 - 2 * p) / (p * (1 - p))
    THAU_2 <- 2 / (p * (1 - p))

    results <- gibbs_sampling(Y, X, initial_values[[5]], initial_values[[2]], THAU_2, THETA, initial_values[[1]], initial_values[[4]], initial_values[[3]], T, N, iterations = ITER, burn = BURN, A_1, A_2, B_1, B_2, C_)

    beta_list <- results[[1]]
    gamma_list <- results[[2]]

    write.csv(beta_list, paste("results/simulations_results/tests/challenge/T", TT, "_t", T, "_beta_p_", p, "_hv_vr_ogz_gc.csv", sep = ""))
    write.csv(gamma_list, paste("results/simulations_results/tests/challenge/T", TT, "_t", T, "_gamma_p_", p, "_hv_vr_ogz_gc.csv", sep = ""))
}

######################## 3
# one date sampler
T <- 180
p <- 0.1
initial_values <- initialize_gibbs(Y, X, T, N, A_1, A_2, B_1, B_2)
THETA <- (1 - 2 * p) / (p * (1 - p))
THAU_2 <- 2 / (p * (1 - p))

results <- gibbs_sampling(Y, X, initial_values[[5]], initial_values[[2]], THAU_2, THETA, initial_values[[1]], initial_values[[4]], initial_values[[3]], T, N, iterations = ITER, burn = BURN, A_1, A_2, B_1, B_2, C_)
################################################################################################################################ 33
smooth_pdf <- function(x, midpoints, pdf_estimates, bandwidth = 0.1) {
    n <- length(midpoints)
    density_est <- 0

    for (i in 1:n) {
        density_est <- density_est + pdf_estimates[i] * dnorm(x, mean = midpoints[i], sd = bandwidth)
    }

    return(density_est)
}
################################################################################################################################ 33
TT <- 50

Y <- as.matrix(read.csv(paste("data/simulations/simulations_T_50_s_5_rho_0.65_sigma_5/y.csv", sep = ""))[, 2])
X <- as.matrix(read.csv(paste("data/simulations/simulations_T_50_s_5_rho_0.65_sigma_5/X.csv", sep = ""))[, 2:31])
x_it <- as.matrix(read.csv(paste("data/simulations/simulations_T_50_s_5_rho_0.65_sigma_5/x_it.csv", sep = ""))[, 2])
Y_it <- as.matrix(read.csv(paste("data/simulations/simulations_T_50_s_5_rho_0.65_sigma_5/y_it.csv", sep = ""))[, 2])
beta_true <- read.csv(paste("data/simulations/simulations_T_50_s_5_rho_0.65_sigma_5/beta_true.csv", sep = ""))[, 2]
N <- ncol(X)
hist(Y_it)
hist(Y)
SFE <- c()
LPS <- c()
p_list <- (1:9) * 0.1

for (T in 45:49) {
    tryCatch(
        {
            y_true <- Y[T + 1]
            beta <- read.csv(paste("results/simulations_results/tests/challenge/T", TT, "_t", T, "_beta_p_0.5_hv_vr_ogz_gc.csv", sep = ""))[, 2:31]
            beta_mean <- colMeans((beta))
            predicted <- X[T + 1, ] %*% beta_mean

            SFE <- append(SFE, (y_true - predicted)^2)
            cat("
        predicted value:", predicted, "--- true value:", y_true)
        },
        error = function(e) {
            # If there's an error (e.g., file not found), print a message
            print(paste("Error reading file:", T, "p: ", p, "- Skipping to next file."))
        }
    )

    quantiles_y_h <- c()
    probs <- c()

    for (p in p_list) {
        tryCatch(
            {
                beta <- read.csv(paste("results/simulations_results/tests/challenge/T", TT, "_t", T, "_beta_p_", p, "_hv_vr_ogz_gc.csv", sep = ""))[, 2:31]
                beta_mean <- apply(beta, 2, mean)

                quantiles_y_h <- append(quantiles_y_h, (beta_mean %*% X[T + 1, ]))
                probs <- append(probs, p)
            },
            error = function(e) {
                # If there's an error (e.g., file not found), print a message
                print(paste("Error reading file:", T, "p: ", p, "- Skipping to next file."))
            }
        )
    }

    quantiles_y_h <- sort(quantiles_y_h)

    delta_q <- diff(quantiles_y_h)
    delta_p <- diff(probs)
    pdf_estimates <- delta_p / delta_q
    midpoints <- (quantiles_y_h[-length(quantiles_y_h)])

    LPS <- append(LPS, log(smooth_pdf(Y[T + 1], midpoints = midpoints, pdf_estimates = pdf_estimates, bandwidth = 2)))
}


MSFE <- mean(SFE)
MSFE
MLPS <- mean(LPS)
MLPS

p <- 0.1
x_it %*% colMeans(read.csv(paste("results/simulations_results/T200_t", T, "_beta_p_", p, ".csv", sep = ""))[, 2:21])
quantile(Y_it, p)

for (T in 40:49) {
    cat("
        ---Period of vintage:", T, "
        ")
    for (p in 1:9) {
        cat("
        ", X[T + 1, ] %*% colMeans(read.csv(paste("results/simulations_results/tests/T", TT, "_t", T, "_beta_p_", p * 0.1, "_u_u_ogz_gc.csv", sep = ""))[, 2:21]))
    }
}
hist(Y_it)

hist(X %*% colMeans(read.csv(paste("results/simulations_results/T200_t", T, "_beta_p_", 0.9, ".csv", sep = ""))[, 2:21]))
set.seed(123)
hist(rep(x_it %*% beta_true, 1000) + rnorm(1000, 0, 1))
hist(Y_it)

################################################################################################################################ 33
TT <- 50

Y <- as.matrix(read.csv(paste("data/simulations/simulations_T_", TT, "/y.csv", sep = ""))[, 2])
X <- as.matrix(read.csv(paste("data/simulations/simulations_T_", TT, "/X.csv", sep = ""))[, 2:21])
x_it <- as.matrix(read.csv(paste("data/simulations/simulations_T_", TT, "/x_it.csv", sep = ""))[, 2])
Y_it <- as.matrix(read.csv(paste("data/simulations/simulations_T_", TT, "/y_it.csv", sep = ""))[, 2])
beta_true <- read.csv(paste("data/simulations/simulations_T_", TT, "/beta_true.csv", sep = ""))[, 2]
N <- ncol(X)
hist(Y_it)
hist(Y)
SFE <- c()
LPS <- c()
probs <- (1:9) * 0.1
beta <- read.csv(paste("results/simulations_results/tests/T", TT, "_t", TT - 10, "_beta_p_0.5_hv_vr_ogz_gc.csv", sep = ""))[, 2:21]
beta_mean <- colMeans((beta))
for (T in 40:49) {
    y_true <- Y[T + 1]
    predicted <- X[T + 1, ] %*% beta_mean

    SFE <- append(SFE, (y_true - predicted)^2)
    cat("
        predicted value:", predicted, "--- true value:", y_true)

    quantiles_y_h <- c()

    for (p in probs) {
        beta <- read.csv(paste("results/simulations_results/tests/T", TT, "_t", TT - 10, "_beta_p_", p, "_hv_vr_ogz_gc.csv", sep = ""))[, 2:21]
        beta_mean <- apply(beta, 2, mean)

        quantiles_y_h <- append(quantiles_y_h, (beta_mean %*% X[T + 1, ]))
    }

    quantiles_y_h <- sort(quantiles_y_h)

    delta_q <- diff(quantiles_y_h)
    delta_p <- diff(probs)
    pdf_estimates <- delta_p / delta_q
    midpoints <- (quantiles_y_h[-length(quantiles_y_h)])

    LPS <- append(LPS, log(smooth_pdf(Y[T + 1], midpoints = midpoints, pdf_estimates = pdf_estimates, bandwidth = 2)))
}

MSFE <- mean(SFE)
MSFE
MLPS <- mean(LPS)
MLPS





colMeans(read.csv(paste("results/simulations_results/tests/T", 50, "/T", 50, "_t", 49, "_beta_p_0.5_u_r_ogz_gc.csv", sep = ""))[, 2:21])