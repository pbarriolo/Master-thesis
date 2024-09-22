#########################################################################
#################### Analysis of results from QR-BMA ####################

# Libraries
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

# Description of data
data_description <- read_excel("C:/Users/pbarr/Documents/ENSAE/3A_MIE/Memoire/forecasting_workspace/data/ALL_DATA.xlsx", sheet = "Data_Description")

##################### Data Import and Stationnarization #################

vd <- "09Q1"
h <- 1 # 0 for nowcasting

df <- get_data(vintage_date = vd, horizon = h)
Y <- as.matrix(na.omit(df[, length(df) - 1])) # predictor training sample
Y_true <- as.matrix(na.omit(df[, length(df)])) # predictor true value vector
X <- (as.matrix(na.omit(df[, 2:(length(df) - 2)]))) # standardize_data() for standardizing
T <- nrow(Y)
apply(X, 2, sd)
sd(Y)
###################### Results Analysis ###############################
path <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/"
setwd(path)

Y_forecast <- read.csv(paste("results/bma/predictions/bma_", vd, "_h", h, "_y_forecast_hv_r_sigma3.csv", sep = ""))
hist(Y_forecast[, 2])
Y_true

########################## Function for smoothing ##########################3
smooth_pdf <- function(x, midpoints, pdf_estimates, bandwidth = 0.1) {
    n <- length(midpoints)
    density_est <- 0

    for (i in 1:n) {
        density_est <- density_est + pdf_estimates[i] * dnorm(x, mean = midpoints[i], sd = bandwidth)
    }

    return(density_est)
}

# compute MSFE on median for QR-BMA
SFE <- c()
LPS <- c()
y_true <- c()
predicted <- c()
h <- 6
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

for (year in c("09", 10:13)) {
    for (quarter in 1:4) {
        if (year >= 15 & quarter >= 1) {
            break()
        }
        vd <- paste(year, "Q", quarter, sep = "")
        df <- get_data(vintage_date = vd, horizon = h)
        Y <- as.matrix(na.omit(df[, length(df) - 1])) # predictor training sample
        Y_true <- as.matrix(na.omit(df[, length(df)])) # predictor true value vector
        X <- as.matrix(na.omit(df[, 2:(length(df) - 2)])) # standardize_data() for standardizing
        T <- nrow(Y)

        beta_mean <- colMeans(read.csv((paste("results/bma/predictions/bma_", vd, "_h", h, "_beta_hv_r_sigma3.csv", sep = "")))[, 2:(ncol(X) + 1)])

        SFE <- append(SFE, (Y_true[length(Y_true)] - X[T + h + 1, ] %*% beta_mean)^2)
        y_true <- append(y_true, Y_true[length(Y_true)])
        predicted <- append(predicted, X[T + h + 1, ] %*% beta_mean)

        # MLPS

        quantiles_y_h <- c()
        Y_forecast <- read.csv(paste("results/bma/predictions/bma_", vd, "_h", h, "_y_forecast_hv_r_sigma3.csv", sep = ""))[, 2]

        for (p in probs) {
            quantiles_y_h <- append(quantiles_y_h, quantile(Y_forecast, p))
        }

        delta_q <- diff(quantiles_y_h)
        delta_p <- diff(probs)
        pdf_estimates <- delta_p / delta_q
        midpoints <- (quantiles_y_h[-length(quantiles_y_h)])

        LPS <- append(LPS, log(smooth_pdf(Y_true[T + 2 * (1 + h)], midpoints = midpoints, pdf_estimates = pdf_estimates, bandwidth = 3)))
    }
}

MSFE <- mean(SFE)
MSFE
MLPS <- mean(na.omit(LPS))
MLPS

plot(1:length(y_true), y_true, type = "l", lty = 1, col = "blue")
lines(1:length(y_true), predicted, type = "l", lty = 1, col = "red")