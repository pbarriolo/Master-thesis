################################################################################
##########    Predictive Densities  from Quantile Regressions   ################
################################################################################

# source: Vulnerable Growth, Adrian, Boyarchenko and Giannone (ABG)

# In this notebook we replicate ABG's parametrics method to estimate probability
# densities from quantile estmates.
# Parameters of a skewed t-distribution are fitted through matching to the estimated
# quantile values.

################################################################################

################### Libraries
library(sn) # for skewed t-distribution functions
library(stats) # for normal distribution
library(pracma) # for non-linear minimization


################### Set path for data
# import data
path <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/results/"
setwd(path)

quantiles_y_h <- read.csv("test3_14Q4_h1_quantile_bma_ns_ns.csv")[2:10]

################## Compute quantiles to match to skewed t-distribution
# forecast =
# q_mean_forecast <- quantiles_y_h[2:10]
# v_q_mean <- c(q_mean_forecast[1], q_mean_forecast[3], q_mean_forecast[5], q_mean_forecast[7], q_mean_forecast[9])

probs <- (1:9) * 0.1
v_probs <- seq(1, 9, 2) * 0.1
v_q <- quantiles_y_h[seq(1, 9, 2)]




# use lsqnonlin to solve the problem

# define funciton to minimize
fn <- function(x) {
  qst(v_probs, xi = x[1], omega = x[2], alpha = x[3], nu = x[4]) - v_q
}
# define initial value
iqn <- qnorm(0.75) - qnorm(0.25)
lc0 <- v_q[1, 5]
sc0 <- 0.5 * (quantiles_y_h[1, 7] + quantiles_y_h[1, 8] - quantiles_y_h[1, 2] - quantiles_y_h[1, 3]) / iqn
sh0 <- 0.
N <- 1.

x0 <- c(lc0, sc0, sh0, N)

result_1 <- lsqnonlin(f, x0)

#### Using optim instead of lsqnonlin

f <- function(x) {
  sum((qst(v_probs, xi = x[1], omega = x[2], alpha = x[3], nu = x[4], tol = 1e-8) - v_q)^2)
}

result <- optim(par = x0, f = f, method = "L-BFGS-B", lower = c(-Inf, 1e-6, -Inf, 0.1), upper = c(Inf, Inf, Inf, Inf))
result


################################################################################
##########    Predictive Densities  from Quantile Regressions   ################
################################################################################

# source: Vulnerable Growth, Adrian, Boyarchenko and Giannone (ABG)

# In this notebook we replicate ABG's parametrics method to estimate probability
# densities from quantile estmates.
# Parameters of a skewed t-distribution are fitted through matching to the estimated
# quantile values.

################################################################################


# libraries
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
path <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/results/"
setwd(path)

vd <- "08Q4"
h <- 1 # 0 for nowcasting

df <- get_data(vintage_date = vd, horizon = h)
Y <- as.matrix(na.omit(df[, length(df) - 1])) # predictor training sample
Y_true <- as.matrix(na.omit(df[, length(df)])) # predictor true value vector
X <- (as.matrix(na.omit(df[, 2:(length(df) - 2)]))) # standardize_data() for standardizing
T <- nrow(Y)
Y_true[length(Y_true)]
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

for (p in c(0.05, (1:9) * 0.1, 0.95)) {
  cat(paste("
    ", p, "-------------"))


  beta <- read.csv(paste("tests/", vd, "_h", h, "_beta_p_", p, "_hv_vr_ogz_gammacorrected.csv", sep = ""))[, 2:19]
  beta_mean <- apply(beta, 2, mean)
  cat(sd(beta[, 1]))
  p_list <- append(p_list, p)

  # gamma_mean <- colMeans(read.csv(paste("predictions/", vd, "_h", h, "_gamma_p_", p * 0.1, ".csv", sep = ""))[, 2:19])
  # cat(gamma_mean)

  quantiles_y_h <- append(quantiles_y_h, (beta_mean %*% X[T + h + 1, ]))
}
# Step 1: Provide Quantiles and Their Probability Levels

# Quantiles (example values)
# quantiles_08 <- c(-8.287388, -1.68651, 1.561752, 3.266263, 5.199962, 7.476578, 9.073098, 10.48642, 14.32321)
# quantiles_14 <- c(-8.33136, -2.108527, 0.1991454, 1.392147, 2.240938, 3.309953, 4.269407, 5.964015, 11.63199)
# quantiles <- c(-4.936641, 0.2723082, 2.055338, 2.633276, 3.235005, 4.594632, 7.434888, 11.2482, 14.0817)
# Corresponding probability levels
probs <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)

# Step 2: Estimate the PDF as the derivative of the CDF

# Compute the differences between quantiles and between probabilities
delta_q <- diff(quantiles_y_h)
delta_p <- diff(probs)

# Estimate the PDF as the derivative of the CDF
pdf_estimates <- delta_p / delta_q

# Associate the midpoints of quantiles with these PDF estimates
midpoints <- (quantiles_y_h[-length(quantiles_y_h)]) # (quantiles[-length(quantiles)] + quantiles[-1])/2

# Step 3: Smooth the Estimated PDF using a Gaussian Kernel

# Define a function to smooth the PDF estimates
smooth_pdf <- function(x, midpoints, pdf_estimates, bandwidth = 0.1) {
  n <- length(midpoints)
  density_est <- 0

  for (i in 1:n) {
    density_est <- density_est + pdf_estimates[i] * dnorm(x, mean = midpoints[i], sd = bandwidth)
  }

  return(density_est)
}

# Step 4: Generate a Range of x Values for Estimating the PDF
x_values <- seq(min(quantiles_y_h) - 1, max(quantiles_y_h) + 1, length.out = 1000)

# Apply the smoothing function to each x value
pdf_values <- sapply(x_values, smooth_pdf, midpoints = midpoints, pdf_estimates = pdf_estimates, bandwidth = 2)

# Step 5: Plot the Smoothed PDF

plot(x_values, pdf_values,
  type = "l", col = "blue", lwd = 2,
  main = "Density one quarter ahead of Inlation in 2008Q4",
  xlab = "x", ylab = "Density"
)
abline(v = Y_true[T + 2 * (1 + h)], col = "red", lty = 2, lwd = 3)

log(smooth_pdf(1, midpoints = midpoints, pdf_estimates = pdf_estimates, bandwidth = 4))


###################################################################################### 3

year_list <- 10:14 # c("08", "09", 10:14)
quarter_list <- 1:4
param_list <- expand.grid(year = year_list, quarter = quarter_list)
h <- 1

LPS <- c()
SFE <- c()
y_true <- c()
predicted <- c()

for (i in 1:nrow(param_list)) {
  # p <- param
  quarter <- param_list[i, "quarter"] #
  year <- param_list[i, "year"] #

  if ((as.numeric(year) < 15) | (quarter == 1)) {
    vd <- paste(year, "Q", quarter, sep = "") #

    df <- get_data(vintage_date = vd, horizon = h) #
    Y <- as.matrix(na.omit(df[, length(df) - 1])) # predictor training sample
    Y_true <- as.matrix(na.omit(df[, length(df)])) # predictor true value vector
    X <- (as.matrix(na.omit(df[, 2:(length(df) - 2)])))
    T <- nrow(Y)

    beta_mean <- colMeans(read.csv((paste("final_predictions/", vd, "_h", h, "_beta_p_0.5_hv_r_ogz_gammacorrected.csv", sep = "")))[seq(1, 5000, 4), 2:(ncol(X) + 1)])
    SFE <- append(SFE, (Y_true[length(Y_true)] - X[T + h + 1, ] %*% beta_mean)^2)
    y_true <- append(y_true, Y_true[length(Y_true)])
    predicted <- append(predicted, X[T + h + 1, ] %*% beta_mean)

    cat("
    predicted value:", X[T + h + 1, ] %*% beta_mean, "--- true value:", Y_true[length(Y_true)])

    quantiles_y_h <- c()
    probs <- c()

    for (p in (1:9) * 0.1) {
      beta <- NULL

      tryCatch(
        {
          beta <- read.csv(paste("final_predictions/", vd, "_h", h, "_beta_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))[, 2:19]
          beta_mean <- apply(beta, 2, mean)

          quantiles_y_h <- append(quantiles_y_h, (beta_mean %*% X[T + h + 1, ]))
          probs <- append(probs, p)
        },
        error = function(e) {
          # If there's an error (e.g., file not found), print a message
          print(paste("Error reading file:", vd, "p: ", p, "- Skipping to next file."))
        }
      )
    }

    quantiles_y_h <- sort(quantiles_y_h)

    delta_q <- diff(quantiles_y_h)
    delta_p <- diff(probs)
    pdf_estimates <- delta_p / delta_q
    midpoints <- (quantiles_y_h[-length(quantiles_y_h)])

    LPS <- append(LPS, log(smooth_pdf(Y_true[T + 2 * (1 + h)], midpoints = midpoints, pdf_estimates = pdf_estimates, bandwidth = 3)))
  } else {
    print("outbound horizon")
  }
}

MLPS <- mean(na.omit(LPS))
MLPS
MSFE <- mean(SFE)
MSFE
plot(1:length(y_true), y_true, type = "l", lty = 1, col = "blue")
lines(1:length(y_true), predicted, type = "l", lty = 1, col = "red")
# read.csv(paste("predictions/10Q1_h1_beta_p_0.5.csv", sep = ""))[, 2:19]