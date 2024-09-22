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

vd <- "08Q4"
h <- 1 # 0 for nowcasting

df <- get_data(vintage_date = vd, horizon = h)
Y <- as.matrix(na.omit(df[, length(df) - 1])) # predictor training sample
Y_true <- as.matrix(na.omit(df[, length(df)])) # predictor true value vector
X <- (as.matrix(na.omit(df[, 2:(length(df) - 2)]))) # standardize_data() for standardizing
T <- nrow(Y)
apply(X, 2, sd)
sd(Y)
###################### Results Analysis ###############################
quantiles_y_h <- c()
p_list <- c()
path <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/results/"
setwd(path)
# p_list <- (1:9) * 0.1
# p_list <- append(p_list, 0.05, 0)
# p_list <- append(p_list, 0.95, 10)
seq(1, 5500, 4)

for (p in (1:9) * 0.1) {
    cat(paste("
    ", p, "-------------"))


    beta <- read.csv(paste("tests/", vd, "_h", h, "_beta_p_", p, "_lv_r_ogz_gammacorrected.csv", sep = ""))[seq(1, 5500, 40), 2:19]
    beta_mean <- apply(beta, 2, mean)
    cat(sd(beta[, 1]))
    p_list <- append(p_list, p)

    # gamma_mean <- colMeans(read.csv(paste("predictions/", vd, "_h", h, "_gamma_p_", p * 0.1, ".csv", sep = ""))[, 2:19])
    # cat(gamma_mean)

    quantiles_y_h <- append(quantiles_y_h, (beta_mean %*% X[T + h + 1, ]))
}

plot(p_list, quantiles_y_h, col = "blue", bg = "blue", pch = 21, main = paste("quantiles for h", h, "and vd", vd, sep = ""))
abline(h = Y_true[T + 2 * (1 + h)], col = "red", lty = 2, lwd = 2)

gamma_ <- read.csv(paste("tests/", vd, "_h", h, "_gamma_p_", 0.5, "_hv_r_ogz_gammacorrected.csv", sep = ""))[, 2:19]
colMeans(gamma_)
plot(1:5500, beta[, 15])

Y_true[T + 2 * (1 + h)] # value we are trying to predict

sd(Y_true)
quantile(Y_true, 0.5)

# beta coefficients posterior visualization
i <- 18 # select the beta we want to observe
for (p in 1:8) {
    beta <- read.csv(paste("tests/test_", vd, "_h", h, "_beta_p_", p * 0.1, "_high_variance_restrictive_ns_s.csv", sep = ""))[, i + 1]
    hist(beta, main = paste("Posterior of Beta", i, "at quantile p=", 0.1 * p), xlab = "Beta")
    cat(median(beta))
}
p <- 5
beta <- read.csv(paste("tests/test_", vd, "_h", h, "_beta_p_", p * 0.1, "_high_variance_restrictive_ns_s.csv", sep = ""))[, i + 1]
hist(beta, main = paste("Posterior of Beta", i, "at quantile p=", 0.1 * p), xlab = "Beta")

# compute MSFE on median for QR-BMA
SFE <- c()
y_true <- c()
predicted <- c()
h <- 1
for (year in 11:14) {
    for (quarter in 1:4) {
        if (year == 15 & quarter > 2) {
            break()
        }
        vd <- paste(year, "Q", quarter, sep = "")
        df <- get_data(vintage_date = vd, horizon = h)
        Y <- as.matrix(na.omit(df[, length(df) - 1])) # predictor training sample
        Y_true <- as.matrix(na.omit(df[, length(df)])) # predictor true value vector
        X <- as.matrix(na.omit(df[, 2:(length(df) - 2)])) # standardize_data() for standardizing
        T <- nrow(Y)

        beta_mean <- colMeans(read.csv((paste("final_predictions/", vd, "_h", h, "_beta_p_0.5_hv_r_ogz_gammacorrected.csv", sep = "")))[seq(1, 5000, 4), 2:(ncol(X) + 1)])

        SFE <- append(SFE, (Y_true[length(Y_true)] - X[T + h + 1, ] %*% beta_mean)^2)
        y_true <- append(y_true, Y_true[length(Y_true)])
        predicted <- append(predicted, X[T + h + 1, ] %*% beta_mean)

        cat("
        predicted value:", X[T + h + 1, ] %*% beta_mean)

        cat("
        --- true value:", y_true)
    }
}

MSFE <- mean(SFE)
MSFE

plot(1:length(y_true), y_true, type = "l", lty = 1, col = "blue")
lines(1:length(y_true), predicted, type = "l", lty = 1, col = "red")

gamma_5 <- read.csv("predictions/15Q1_h1_gamma_p_0.5_u_u_v2.csv")[1000:11000, 2:19]
gamma_52 <- read.csv("predictions/15Q1_h1_gamma_p_0.5_hv_r_v2.csv")[1000:11000, 2:19]
colMeans(gamma_5)
colMeans(gamma_52)