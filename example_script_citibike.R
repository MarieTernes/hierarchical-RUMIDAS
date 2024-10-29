###################################################################################################################################
#------- Example script for obtaining the hierarchical estimator for pooled RUMIDAS with the citibike application -----------------#
# "Hierarchical Regularizers for Reverse Unrestricted Mixed Data Sampling Regressions" by Alain Hecq, Marie Ternes and Ines Wilms  #
###################################################################################################################################
rm(list=ls())     # Clean memory

# Load libraries
library(zoo)
library(RSpectra)

# Set working directory to the folder where this script is contained in
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Main functions
source("RUMIDAS_functions.R")
source("sparseRUMIDAS_functions.R")

#########################################################################
#-------------- Data example & Data pre-processing ---------------------#
#########################################################################
# Load Data
setwd("./data")
load("citibike_rumidas.RData")
load("rides.RData")
names(rides)
# [1] "Date"                "Subways"             "Buses"               "LIRR"                "Metro"              
# [6] "AccessARide"         "BridgesAndTunnes"    "StatenIslandRailway"

# Data series
area = "NYC"
x_hf = citibike_rumidas$m2[,area] # hourly data (high-frequency)
area = "StatenIslandRailway"
get_col = which(names(rides)==area)
y_lf = rides[, get_col] # daily data (low-frequency)
datetime <- as.POSIXct(paste(citibike_rumidas$m2$year,citibike_rumidas$m2$month,citibike_rumidas$m2$day,citibike_rumidas$m2$hour, sep = "-"), format = "%Y-%m-%d-%H", tz = "America/New_York")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#####################################
#----------- Parameters -------------#
######################################
m = 24 #frequency mismatch
p_hf = m*7 # number of HF lags
p_lf = 7 # number of LF lags
l_lambda = 20 # number of lambdas in grid
h = 1 # horizon

window_days = 30 #
tin = 1:(window_days*m)
tout = window_days*m+h



###########################################################
#-------------- Pooled RUMIDAS model ---------------------#
###########################################################

######################################
#--------- Data preparation ---------#
######################################
#### Prepare data for (mixed-frequency) lags
## IF YOU USE THE POOLED RUMIDAS USE flexLF = T, FOR THE ORIGINAL RUMIDAS flexLF = F
algn = align_data(x_hf = x_hf, y_lf = y_lf, p_hf = p_hf, p_lf = p_lf, m = m, h = h, flexLF = T)

# other input parameters
inp_mflags = inputs_algorithm(algn$k_hf, algn$k_lf, p_hf, p_lf, m, penalty = "HIER", flexLF = T)
group_index_mflags = inp_mflags$group_index
penalty_mflags = inp_mflags$penalty
columns_LF_mflags = inp_mflags$columns_LF
columns_HF_mflags = inp_mflags$columns_HF
# the inp_mflags$penalty_values returns are the recency-based penalty values: 
# penalty_values_mflags = inp_mflags$penalty_values
# however it is also possble to hard-code the penalty values if one has a specific penalty structure in mind
# below we use a seasonality based penalty structure 

lfpen_element671 = c(2:p_lf, 1)
lfpen_element617 = c(2:(p_lf-1), 1, p_lf)
# HARD-CODED LINE:
if(h == 1){
  LFpen = rep(c(rep(c(c(2:p_lf), 1), m)), each = algn$k_lf)
  HFpen = rep(c(2:(p_hf), 1), each = algn$k_hf)    
  penalty_values_mflags = c(LFpen, HFpen) 
}
if(h>1 & h <= 24){
  LFpen = rep( c(rep(lfpen_element671, m-h+1), rep(lfpen_element617, h-1)) , each = algn$k_lf)
  HFpen = rep(c(2:(p_hf-h+1),1, (p_hf-h+2):p_hf), each = algn$k_hf)
  penalty_values_mflags = c(LFpen, HFpen) 
}
if(h == 24*7){
  LFpen = rep(c(c(2:p_lf,1), c(rep(c(c(2:(p_lf-1)), 1, p_lf), m-1))), each = algn$k_lf)
  HFpen = rep(c(1:p_hf), each = algn$k_hf)    
  penalty_values_mflags = c(LFpen, HFpen) 
}


#####------------------------------ Hierarchical Lasso ------------------------------------######

# Split into in and out-of-sample

ymflags_in = algn$y[tin]
Xmflags_in = algn$X[tin,]
ymflags_out = algn$y[tout]
Xmflags_out = algn$X[tout,]

##### before penalization still standardize!! #####
mflags_s = rumidas_standardize_flexLF(y = ymflags_in, Xdummy_LF = Xmflags_in[,columns_LF_mflags], Xmflags_in[,columns_HF_mflags], m)
Xmflags_ins = mflags_s$Xd_s
ymflags_ins = mflags_s$yd_s

mu_mflags_y = mflags_s$mu_y
sigma_mflags_y = mflags_s$sigma_y
mu_mflags_X = mflags_s$mu_X
sigma_mflags_X = mflags_s$sigma_X


#############################################################################
#-------------- Obtaining the hierarchical regularizer ---------------------#
#############################################################################

########## POOLED RUMIDAS ###############
######## TUNING PARAMETER SELECTION BY BIC ##########
# mflags_hier: hierarchical on LF and HF
mflags_hier_out = sparseRUMIDAS(y = ymflags_ins, X = Xmflags_ins, group_index = group_index_mflags, penalty = penalty_mflags, penalty_values = penalty_values_mflags, 
                                  columns_LF = columns_LF_mflags, columns_HF = columns_HF_mflags,
                                  l_lambda = l_lambda, post_lasso = TRUE, epsilon = 10^-4, max_iter = 500)
betas_mflags_hier = mflags_hier_out$betas[,mflags_hier_out$bic_index]
lambda_bic = mflags_hier_out$lambdaSeq[mflags_hier_out$bic_index]


## Forecast (re-transform back after standardization)
# standardize X for HIER out-of-sample
index_forecast_mflags = which(Xmflags_out!=0)
Xmflags_outs = ( Xmflags_out[index_forecast_mflags] - mu_mflags_X[index_forecast_mflags] ) / sigma_mflags_X[index_forecast_mflags]
ymflags_hier_hats = Xmflags_outs %*% betas_mflags_hier[index_forecast_mflags]
ymflags_hier_hat = sigma_mflags_y * ymflags_hier_hats + mu_mflags_y

# error
ymflags_out - max(ymflags_hier_hat,0)


######## ALTERNATIVE: TUNING PARAMETER SELECTION BY TIME SERIES CROSS-VALIDATION ##########
# Warning takes long 
mflags_hier_out = sparseRUMIDAS_tscv(y = ymflags_ins, X = Xmflags_ins, cvcut = 168, 
                                     group_index = group_index_mflags, penalty = penalty_mflags, penalty_values = penalty_values_mflags, 
                                     columns_LF = columns_LF_mflags, columns_HF = columns_HF_mflags,
                                     l_lambda = l_lambda, post_lasso = TRUE, epsilon = 10^-4, max_iter = 500)
betas_mflags_hier = mflags_hier_out$betas
lambda_tscv_min = mflags_hier_out$lambda_min
CVscore_save = mflags_hier_out$MAFE_avg

