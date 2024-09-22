# Quantile Regression Bayesian Model Averaging
# author: Pablo Barrio

# libraries
library(readxl)
library(zoo)
library(tidyverse)
library(mvtnorm)
library(ghyp)
library(robustbase)

# sources for function
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_transform.R")
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_gibbs.R")
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_import_data.R")
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_ridge.R")


# Import df into one df_frame


# Description of data
data_description <- read_excel("C:/Users/pbarr/Documents/ENSAE/3A_MIE/Memoire/forecasting_workspace/data/ALL_DATA.xlsx", sheet = "Data_Description")

##################### Data Import and Stationnarization #################

vd <- "14Q4"
h <- 6 # 0 for nowcasting

df <- get_data(vintage_date = vd, horizon = h)
Y <- as.matrix(na.omit(df[, length(df) - 1])) # predictor training sample
Y_true <- as.matrix(na.omit(df[, length(df)])) # predictor true value vector
X <- as.matrix(na.omit(df[, 2:(length(df) - 2)])) # standardize_data() for standardizing


############################## Initialization  ###################################
##############   code for only one vintage data and horizons  ####################

# Variables and Parameters definition

ITER <- 200
BURN <- 100

N <- ncol(X) # number of covariables
C_ <- 0.00001 # try with 0
A_1 <- 2 # test with high variance (1, 10) and medium variance (2, 1), (2, 0.1) high variance, uninformative (10,10)
A_2 <- 0.5
B_1 <- 1 # ?? test with (5,1) loose, (1,5) restrictive, and (1, 20) very restrictive, uninformative (1,1)
B_2 <- 1 # ??

T <- nrow(Y) # training sample size, depends on the vintage and the forecast horizon

path <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/results/all_sample/"

# Initiate variables, parameters and hyperparameters

initial_values <- initialize_gibbs(Y, X, T, N, B_1, B_2)

################### 3######## Run Gibbs ###################################

for (p in (1:9) * 0.1) {
  THETA <- (1 - 2 * p) / (p * (1 - p))
  THAU_2 <- 2 / (p * (1 - p))

  beta_list <- gibbs_sampling(Y, X, initial_values[[5]], initial_values[[2]], THAU_2, THETA, initial_values[[1]], initial_values[[4]], initial_values[[3]], T, N, iterations = ITER, burn = BURN, A_1, A_2, B_1, B_2, C_, h, path, p, vd, keep_track = TRUE)
}
# beta_list_test <- gibbs_sampling_prediction(Y, X, initial_values[[5]], initial_values[[2]], THAU_2, THETA, initial_values[[1]], initial_values[[4]], initial_values[[3]], T, N, iterations = ITER, burn = BURN, A_1, A_2, B_1, B_2, C_, h, path, p, vd, keep_track = TRUE)


#################### Predictions for 2010Q1 to 2015q4 for horizons 1 and 6
year_list <- c(14)
quarter_list <- 2:4
path_2 <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/"
setwd(path_2)

ITER <- 5000
BURN <- 500
for (year in year_list) {
  for (quarter in quarter_list) {
    for (h in c(1, 6)) {
      if ((h == 1) | (h == 6 & (year < 15 | quarter <= 1))) {
        vintage_date <- paste(year, "Q", quarter, sep = "")
        subset_raw_data <- select_data(raw_data, tcode, vintage_date, h)
        # transform series
        df <- transform_data(subset_raw_data[2:nrow(subset_raw_data), ], tcode)
        #  predictor training sample
        Y <- as.matrix(na.omit(df[, length(df) - 1]))
        # covariables training sample
        X <- as.matrix(na.omit(df[, 2:(length(df) - 2)]))
        # predictor true value vector
        Y_true <- as.matrix(na.omit(df[, length(df)]))

        T <- nrow(Y) # training sample size, depends on the vintage and the forecast horizon

        mu <- 0.5

        # Initiate variables, parameters and hyperparameters

        beta_0 <- rep(0, ncol(X)) # ridge_regression(Y, X, mu, N) # should I be using a quantile regression??

        var_residuals <- t(Y - X[1:T, ] %*% beta_0) %*% (Y - X[1:T, ] %*% beta_0) / (T - N - 1)

        delta_0 <- diag(var_residuals[1] * solve(t(X[1:T, ]) %*% X[1:T, ] + mu * diag(N)) %*% t(X[1:T, ]) %*% X[1:T, ] %*% solve(t(X[1:T, ]) %*% X[1:T, ] + mu * diag(N)))

        gamma_0 <- initiate_gamma(beta_0)

        Z_0 <- rep(1, T) # rexp(T)*0

        for (p in 1:9 * 0.1) {
          cat("
              Gibbs Sampler for p=", p, "for date:", vintage_date, "and horizon:", h)

          THETA <- (1 - 2 * p) / (p * (1 - p))
          THAU_2 <- 2 / (p * (1 - p))
          cat("
                  theta:", THETA, "
                  thau_2:", THAU_2)

          beta_list <- gibbs_sampling_prediction(Y, X, Z_0, pi_0, THAU_2, THETA, beta_0, gamma_0, delta_0, T, N, iterations = ITER, burn = BURN, A_1, A_2, B_1, B_2, C_, h, path, p, vintage_date, keep_track = TRUE)

          write.csv(beta_list, paste("results/", vintage_date, "_h", h, "_beta_p_", p, ".csv", sep = ""))
        }
      }
    }
  }
}

year <- 11
quarter <- 4
h <- 6

vintage_date <- paste(year, "Q", quarter, sep = "")
subset_raw_data <- select_data(raw_data, tcode, vintage_date, h)
# transform series
df <- transform_data(subset_raw_data[2:nrow(subset_raw_data), ], tcode)
#  predictor training sample
Y <- as.matrix(na.omit(df[, length(df) - 1]))
# covariables training sample
X <- as.matrix(na.omit(df[, 2:(length(df) - 2)]))
# predictor true value vector
Y_true <- as.matrix(na.omit(df[, length(df)]))

T <- nrow(Y) # training sample size, depends on the vintage and the forecast horizon

mu <- 0.5

# Initiate variables, parameters and hyperparameters

beta_0 <- rep(0, ncol(X)) # ridge_regression(Y, X, mu, N) # should I be using a quantile regression??

var_residuals <- t(Y - X[1:T, ] %*% beta_0) %*% (Y - X[1:T, ] %*% beta_0) / (T - N - 1)

delta_0 <- diag(var_residuals[1] * solve(t(X[1:T, ]) %*% X[1:T, ] + mu * diag(N)) %*% t(X[1:T, ]) %*% X[1:T, ] %*% solve(t(X[1:T, ]) %*% X[1:T, ] + mu * diag(N)))

gamma_0 <- initiate_gamma(beta_0)

Z_0 <- rep(1, T) # rexp(T)*0

for (p in 1:9 * 0.1) {
  cat("
              Gibbs Sampler for p=", p, "for date:", vintage_date, "and horizon:", h)

  THETA <- (1 - 2 * p) / (p * (1 - p))
  THAU_2 <- 2 / (p * (1 - p))
  cat("
                  theta:", THETA, "
                  thau_2:", THAU_2)

  beta_list <- gibbs_sampling_prediction(Y, X, Z_0, pi_0, THAU_2, THETA, beta_0, gamma_0, delta_0, T, N, iterations = ITER, burn = BURN, A_1, A_2, B_1, B_2, C_, h, path, p, vintage_date, keep_track = TRUE)

  write.csv(beta_list, paste("results/", vintage_date, "_h", h, "_beta_p_", p, ".csv", sep = ""))
}