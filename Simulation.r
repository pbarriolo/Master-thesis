####### Data simulation #######
###############################

# We simulate data to study properties of QR-BMA
# author: Pablo Barrio

# libraries
library(readxl)
library(zoo)
library(tidyverse)
library(mvtnorm)
library(ghyp)
library(MASS)
library(robustbase)

# sources for function
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_transform.R")
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_gibbs.R")
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_import_data.R")
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_ridge.R")
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_simulations.R")

# Import df into one df_frame
path <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/data/simulations/"
setwd(path)

T_list <- c(200, 100, 50)

# We first simulate data with different sample sizes
# other parameters are not moved
for (T in T_list) {
    # set up other parameters
    k <- 20
    s <- 10
    rho <- 0.75

    # run simulation
    results_sim <- simulate_data_one_value(T, k, s, rho, 1)
    y <- results_sim[[1]]
    X <- results_sim[[2]]
    beta_true <- results_sim[[3]]
    y_it <- results_sim[[4]]
    x_it <- results_sim[[5]]

    write.csv(y, paste("simulations_T_", T, "/y.csv", sep = ""))
    write.csv(X, paste("simulations_T_", T, "/X.csv", sep = ""))
    write.csv(beta_true, paste("simulations_T_", T, "/beta_true.csv", sep = ""))
    write.csv(y_it, paste("simulations_T_", T, "/y_it.csv", sep = ""))
    write.csv(x_it, paste("simulations_T_", T, "/x_it.csv", sep = ""))
}

# we run simulations for a hard dataset that wil challenge the model
hist(y_it)
T <- 50 # medium size for sample
k <- 30
s <- 5 # low number of true predictors
rho <- 0.65
SIGMA <- 5 # high variance of the error term to make prediction more challenging

results_sim <- simulate_data_one_value(T, k, s, rho, SIGMA, iterations = 1000)
y <- results_sim[[1]]
X <- results_sim[[2]]
beta_true <- results_sim[[3]]
y_it <- results_sim[[4]]
x_it <- results_sim[[5]]

write.csv(y, paste("simulations_T_50_s_5_rho_0.65_sigma_5/y.csv", sep = ""))
write.csv(X, paste("simulations_T_50_s_5_rho_0.65_sigma_5/X.csv", sep = ""))
write.csv(beta_true, paste("simulations_T_50_s_5_rho_0.65_sigma_5/beta_true.csv", sep = ""))
write.csv(y_it, paste("simulations_T_50_s_5_rho_0.65_sigma_5/y_it.csv", sep = ""))
write.csv(x_it, paste("simulations_T_50_s_5_rho_0.65_sigma_5/x_it.csv", sep = ""))