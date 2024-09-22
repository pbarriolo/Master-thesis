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

# Import df into one df_frame
path <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/"
setwd(path)

# set up cores for parallel computation
numCores <- detectCores()
registerDoParallel(3)

# Description of data
data_description <- read_excel("C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/data/ALL_DATA.xlsx", sheet = "Data_Description")

##################### Data Import and Stationnarization #################
year <- "09"
quarter <- 4
vd <- paste(year, "Q", quarter, sep = "")
h <- 1 # 0 for nowcasting

df <- get_data(vintage_date = vd, horizon = h)
Y <- as.matrix(na.omit(df[, length(df) - 1])) #  predictor training sample
Y_true <- (as.matrix(na.omit(df[, length(df)]))) # predictor true value vector
X <- as.matrix(na.omit(df[, 2:(length(df) - 2)])) # standardize_data() for standardizing


############################## Initialization  ###################################
##############   code for only one vintage data and horizons  ####################

# Variables and Parameters definition
ITER <- 5000
BURN <- 500

N <- ncol(X) # number of covariables
T <- nrow(Y) # training sample size, depends on the vintage and the forecast horizon
C_ <- 0.00001 # try with 0
A_1 <- 1 # high variance hv  (1, 10) and medium variance mv (2, 1), (2, 0.1) low variance lv , uninformative (0.1,0.1) u
A_2 <- 10
B_1 <- 1 # (5,1) loose l , (1,5) restrictive r, and (1, 20) very restrictive vr, uninformative (1,1) u
B_2 <- 5
mu <- 0.5

# initial_values <- initialize_gibbs(Y, X, T, N, A_1, A_2, B_1, B_2)
year_list <- c("08", "09", 10:14)
p_list <- (1:9) * 0.1
# p_list <- c(0.05, 0.95)
# p_list <- append(p_list, c(0.05, 0.95))
quarter_list <- c(1:4)
h_list <- c(6)
param_list <- expand.grid(year = year_list, quarter = quarter_list, p = p_list, h = h_list)
# Import df into one df_frame
path <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/"
setwd(path)

foreach(i = 1:nrow(param_list), .combine = c, .packages = c("readxl", "zoo", "tidyverse", "mvtnorm", "ghyp", "robustbase", "dplyr")) %dopar% {
    # p <- param
    p <- param_list[i, "p"] #
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

    initial_values <- initialize_gibbs(Y, X, T, N, A_1, A_2, B_1, B_2) #

    if ((h == 1 & (as.numeric(year) < 15 | quarter <= 2)) | (h == 6 & (as.numeric(year) < 14 | (as.numeric(year) == 14 & quarter <= 1)))) {
        THETA <- (1 - 2 * p) / (p * (1 - p))
        THAU_2 <- 2 / (p * (1 - p))

        results <- gibbs_sampling(Y, X, initial_values[[5]], initial_values[[2]], THAU_2, THETA, initial_values[[1]], initial_values[[4]], initial_values[[3]], T, N, iterations = ITER, burn = BURN, A_1, A_2, B_1, B_2, C_)

        beta_list <- results[[1]]
        gamma_list <- results[[2]]
        delta_list <- results[[3]]
        Z_list <- results[[4]]
        pi_list <- results[[5]]
        l_0_list <- results[[6]]
        l_1_list <- results[[7]]

        write.csv(beta_list, paste("results/final_predictions/", vd, "_h", h, "_beta_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
        # write.csv(gamma_list, paste("results/final_predictions/", vd, "_h", h, "_gamma_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
        # write.csv(delta_list, paste("results/final_predictions/", vd, "_h", h, "_delta_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
        # write.csv(Z_list, paste("results/final_predictions/", vd, "_h", h, "_Z_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
        # write.csv(pi_list, paste("results/final_predictions/", vd, "_h", h, "_pi0_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
        # write.csv(l_0_list, paste("results/final_predictions/", vd, "_h", h, "_l0_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
        # write.csv(l_1_list, paste("results/final_predictions/", vd, "_h", h, "_l1_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
    }
}

p <- 0.5
beta_list <- read.csv(paste("results/all_sample/", vd, "_h", h, "_beta_p_", p, "_u_l_ogz_gammacorrected.csv", sep = ""))[, 2:19]
gamma_list <- read.csv(paste("results/all_sample/", vd, "_h", h, "_gamma_p_", p, "_u_l_ogz_gammacorrected.csv", sep = ""))[, 2:19]
delta_list <- read.csv(paste("results/all_sample/", vd, "_h", h, "_delta_p_", p, "_u_l_ogz_gammacorrected.csv", sep = ""))[, 2:19]
Z_list <- read.csv(paste("results/all_sample/", vd, "_h", h, "_Z_p_", p, "_u_l_ogz_gammacorrected.csv", sep = ""))
pi_0_list <- read.csv(paste("results/all_sample/", vd, "_h", h, "_pi0_p_", p, "_u_l_ogz_gammacorrected.csv", sep = ""))[, 2]
l_0_list <- read.csv(paste("results/all_sample/", vd, "_h", h, "_l0_p_", p, "_u_l_ogz_gammacorrected.csv", sep = ""))[, 2]
l_1_list <- read.csv(paste("results/all_sample/", vd, "_h", h, "_l1_p_", p, "_u_l_ogz_gammacorrected.csv", sep = ""))[, 2]


# when need to recompute for one quantile
year <- "10"
quarter <- 2
vd <- paste(year, "Q", quarter, sep = "")
h <- 1 # 0 for nowcasting
df <- get_data(vintage_date = vd, horizon = h)
Y <- (as.matrix(na.omit(df[, length(df) - 1]))) # predictor training sample
Y_true <- as.matrix(na.omit(df[, length(df)])) # predictor true value vector
X <- (as.matrix(na.omit(df[, 2:(length(df) - 2)]))) # standardize_data() for standardizing
ITER <- 5000
BURN <- 500
N <- ncol(X) # number of covariables
T <- nrow(Y) # training sample size, depends on the vintage and the forecast horizon
C_ <- 0.00001 # try with 0
A_1 <- 1 # high variance hv  (1, 10) and medium variance mv (2, 1), (2, 0.1) low variance lv , uninformative (0.1,0.1) u
A_2 <- 10
B_1 <- 1 # (5,1) loose l , (1,5) restrictive r, and (1, 20) very restrictive vr, uninformative (1,1) u
B_2 <- 5
mu <- 0.5
initial_values <- initialize_gibbs(Y, X, T, N, A_1, A_2, B_1, B_2)
p <- 0.5
THETA <- (1 - 2 * p) / (p * (1 - p))
THAU_2 <- 2 / (p * (1 - p))
#
results <- gibbs_sampling(Y, X, initial_values[[5]], initial_values[[2]], THAU_2, THETA, initial_values[[1]], initial_values[[4]], initial_values[[3]], T, N, iterations = ITER, burn = BURN, A_1, A_2, B_1, B_2, C_)

beta_list <- results[[1]]
gamma_list <- results[[2]]
delta_list <- results[[3]]
Z_list <- results[[4]]
pi_0_list <- results[[5]]
l_0_list <- results[[6]]
l_1_list <- results[[7]]

sequ <- seq(1, 110000, 8)
par(mfrow = c(2, 1))
plot(1:5500, log(l_0_list[, 1]), type = "l", main = "Trace of li_0")
plot(1:5500, log(l_1_list[, 1]), type = "l", main = "Trace of li_1")

el <- 1
par(mfrow = c(2, 2))
plot(beta_list[, el], main = "Trace of Beta[1]")
plot(gamma_list[, el], main = "Trace of Gamma[1]")
plot(Z_list[, T], type = "l", main = "Trace of Z[1]")
plot(delta_list[, el], type = "l", main = "Trace of delta")
colMeans(gamma_list)

X[T + 1 + h, ] %*% colMeans(beta_list)

write.csv(beta_list, paste("results/final_predictions/", vd, "_h", h, "_beta_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
# write.csv(gamma_list, paste("results/final_predictions/", vd, "_h", h, "_gamma_p_", p, "_hv_r_og.csv", sep = ""))
# write.csv(delta_list, paste("results/final_predictions/", vd, "_h", h, "_delta_p_", p, "_hv_r_og.csv", sep = ""))
# write.csv(Z_list, paste("results/final_predictions/", vd, "_h", h, "_Z_p_", p, "_hv_r_og.csv", sep = ""))
# write.csv(pi_list, paste("results/final_predictions/", vd, "_h", h, "_pi0_p_", p, "_hv_r_og.csv", sep = ""))
# write.csv(l_0_list, paste("results/final_predictions/", vd, "_h", h, "_l0_p_", p, "_hv_r_og.csv", sep = ""))
# write.csv(l_1_list, paste("results/final_predictions/", vd, "_h", h, "_l1_p_", p, "_hv_r_og.csv", sep = ""))

# commentaire: il semblerait qu'on a plus les memes resultats d'avant.
# la variance de la prediction, difference entre les quantiles a diminué
# je sais pas si bonne nouvelle ou mauvaise
# en tout cas pb, car pas encore capable de predire le bon chiffre
# median est mal predite
# possibles raisons du changement de resultats: changement de posterior de z ou changement dqns initialisation
# très probablement à cause du changement de posterior de Z

ITER <- 5000
BURN <- 500
C_ <- 0.00001 # try with 0
A_1 <- 1 # high variance hv  (1, 10) and medium variance mv (2, 1), (2, 0.1) low variance lv , uninformative (0.1,0.1) u
A_2 <- 10
B_1 <- 1 # (5,1) loose l , (1,5) restrictive r, and (1, 20) very restrictive vr, uninformative (1,1) u
B_2 <- 5

# initial_values <- initialize_gibbs(Y, X, T, N, A_1, A_2, B_1, B_2)
year_list <- c("08") # , "09") # , 10:14)
p_list <- (1:9) * 0.1
# p_list <- c(0.05, 0.95)
# p_list <- append(p_list, c(0.05, 0.95))
quarter_list <- c(3)
h_list <- c(1)
param_list <- expand.grid(year = year_list, quarter = quarter_list, p = p_list, h = h_list)
# Import df into one df_frame
path <- "C:/Users/pbarr/OneDrive/Estudios/ENSAE/3A_MIE/Memoire/forecasting_workspace/"
setwd(path)

for (i in 1:nrow(param_list)) {
    # p <- param
    p <- param_list[i, "p"] #
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

    initial_values <- initialize_gibbs(Y, X, T, N, A_1, A_2, B_1, B_2) #

    if ((h == 1 & (as.numeric(year) < 15 | quarter <= 2)) | (h == 6 & (as.numeric(year) < 14 | (as.numeric(year) == 14 & quarter <= 1)))) {
        THETA <- (1 - 2 * p) / (p * (1 - p))
        THAU_2 <- 2 / (p * (1 - p))

        results <- gibbs_sampling(Y, X, initial_values[[5]], initial_values[[2]], THAU_2, THETA, initial_values[[1]], initial_values[[4]], initial_values[[3]], T, N, iterations = ITER, burn = BURN, A_1, A_2, B_1, B_2, C_)

        beta_list <- results[[1]]
        gamma_list <- results[[2]]
        delta_list <- results[[3]]
        Z_list <- results[[4]]
        pi_list <- results[[5]]
        l_0_list <- results[[6]]
        l_1_list <- results[[7]]

        write.csv(beta_list, paste("results/final_predictions/", vd, "_h", h, "_beta_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
        # write.csv(gamma_list, paste("results/final_predictions/", vd, "_h", h, "_gamma_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
        # write.csv(delta_list, paste("results/final_predictions/", vd, "_h", h, "_delta_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
        # write.csv(Z_list, paste("results/final_predictions/", vd, "_h", h, "_Z_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
        # write.csv(pi_list, paste("results/final_predictions/", vd, "_h", h, "_pi0_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
        # write.csv(l_0_list, paste("results/final_predictions/", vd, "_h", h, "_l0_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
        # write.csv(l_1_list, paste("results/final_predictions/", vd, "_h", h, "_l1_p_", p, "_hv_r_ogz_gammacorrected.csv", sep = ""))
    }
}
as.numeric("08")