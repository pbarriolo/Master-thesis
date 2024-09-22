# Paralelisation of gibbs sampler
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
source("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/functions/functions_gibbs_BMA.R")

# Import df into one df_frame
path <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/"
setwd(path)

# set up cores for parallel computation
numCores <- detectCores()
registerDoParallel(3)

# Description of data
data_description <- read_excel("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/data/ALL_DATA.xlsx", sheet = "Data_Description")

##################### Data Import and Stationnarization #################
year <- "14"
quarter <- 1
vd <- paste(year, "Q", quarter, sep = "")
h <- 6 # 0 for nowcasting

df <- get_data(vintage_date = vd, horizon = h)
Y <- as.matrix(na.omit(df[, length(df) - 1])) # predictor training sample
Y_true <- as.matrix(na.omit(df[, length(df)])) # predictor true value vector
X <- as.matrix(na.omit(df[, 2:(length(df) - 2)])) # standardize_data() for standardizing


############################## Initialization  ###################################
##############   code for only one vintage data and horizons  ####################

# Variables and Parameters definition
ITER <- 5000
BURN <- 500

N <- ncol(X) # number of covariables
T <- nrow(Y) # training sample size, depends on the vintage and the forecast horizon
C_ <- 0.00001 # try with 0
A_1 <- 1 # high variance hv  (1, 10) and medium variance mv (2, 1), (2, 0.1) low variance lv , uninformative (10,10) u
A_2 <- 10
B_1 <- 1 # (5,1) loose l , (1,5) restrictive r, and (1, 20) very restrictive vr, uninformative (1,1) u
B_2 <- 5
SIGMA_2 <- 3

# benchmark has uninformative priors

pi_0 <- 0.1

initial_values <- initialize_gibbs_bma(Y, X, T, N, B_1, B_2)
year_list <- 10:15
quarter_list <- 1:4
param_list <- expand.grid(year = year_list, quarter = quarter_list)


results <- gibbs_sampling_bma(Y, X, SIGMA_2, initial_values[[1]], initial_values[[2]], initial_values[[3]], pi_0, T, N, iterations = ITER, burn = BURN, A_1, A_2, B_1, B_2, C_)

beta_list <- results[[1]]
gamma_list <- results[[2]]
Y_forecast <- results[[3]]

write.csv(beta_list, paste("results/bma/predictions/bma_", vd, "_h", h, "_beta_u_r_sigma1.csv", sep = ""))
write.csv(gamma_list, paste("results/bma/predictions/bma_", vd, "_h", h, "_gamma_u_r_sigma1.csv", sep = ""))
write.csv(Y_forecast, paste("results/bma/predictions/bma_", vd, "_h", h, "_y_forecast_u_r_sigma1.csv", sep = ""))

hist(Y_forecast[, 1])
Y_true


############################################################################################

# Variables and Parameters definition
ITER <- 5000
BURN <- 500
C_ <- 0.00001 # try with 0
A_1 <- 1 # high variance hv  (1, 10) and medium variance mv (2, 1), (2, 0.1) low variance lv , uninformative (10,10) u
A_2 <- 10
B_1 <- 1 # (5,1) loose l , (1,5) restrictive r, and (1, 20) very restrictive vr, uninformative (1,1) u
B_2 <- 5
SIGMA_2 <- 3
# initial_values <- initialize_gibbs(Y, X, T, N, A_1, A_2, B_1, B_2)
year_list <- c("08", "09", 10:15)
quarter_list <- 1:4
h_list <- c(1, 6)
param_list <- expand.grid(year = year_list, quarter = quarter_list, h = h_list)
# Import df into one df_frame
path <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/"
setwd(path)

foreach(i = 1:nrow(param_list), .combine = c, .packages = c("readxl", "zoo", "tidyverse", "mvtnorm", "ghyp", "robustbase", "dplyr")) %dopar% {
    # p <- param
    quarter <- param_list[i, "quarter"] #
    year <- param_list[i, "year"] #
    h <- param_list[i, "h"]

    vd <- paste(year, "Q", quarter, sep = "") #

    df <- get_data(vintage_date = vd, horizon = h) #
    Y <- as.matrix(na.omit(df[, length(df) - 1])) # predictor training sample
    Y_true <- as.matrix(na.omit(df[, length(df)])) # predictor true value vector
    X <- as.matrix(na.omit(df[, 2:(length(df) - 2)])) # standardize_data() for standardizing
    T <- nrow(Y) #
    N <- ncol(X) #

    initial_values <- initialize_gibbs_bma(Y, X, T, N, B_1, B_2)

    if ((h == 1 & (as.numeric(year) < 15 | quarter <= 2)) | (h == 6 & (as.numeric(year) < 14 | (as.numeric(year) == 14 & quarter <= 1)))) {
        results <- gibbs_sampling_bma(Y, X, SIGMA_2, initial_values[[1]], initial_values[[2]], initial_values[[3]], pi_0, T, N, iterations = ITER, burn = BURN, A_1, A_2, B_1, B_2, C_)

        beta_list <- results[[1]]
        # gamma_list <- results[[2]]
        Y_forecast <- results[[3]]

        write.csv(beta_list, paste("results/bma/predictions/bma_", vd, "_h", h, "_beta_hv_r_sigma3.csv", sep = ""))
        # write.csv(gamma_list, paste("results/bma/predictions/bma_", vd, "_h", h, "_gamma_u_r_sigma1.csv", sep = ""))
        write.csv(Y_forecast, paste("results/bma/predictions/bma_", vd, "_h", h, "_y_forecast_hv_r_sigma3.csv", sep = ""))
    }
}