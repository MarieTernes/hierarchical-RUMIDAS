###################################################################################################################################
#--- Example script for obtaining the hierarchical estimator for pooled RUMIDAS and RUMIDAS with the volatility application -------#
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
load("SP500_approx20.RData")
load("FREDMDtransxf20.RData")

# Data series
LF_select = data.frame(INDPRO = FREDMDtransxf20$INDPRO)
x_hf = log(SP500_approx20$rv5*10000) #log(SP500_approx20$medrv*10000) 
y_lf = as.matrix(LF_select)
macrovar <- "indpro"

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#####################################
#----------- Parameters -------------#
######################################
m = 20 #frequency mismatch
p_hf = m # number of HF lags
p_lf = 1 # number of LF lags
l_lambda = 20 # number of lambdas in grid
h = 1 # horizon

# window_days = 30 # 5*28 # 30 # 5*28
# tin = 1:(window_days*m)
# tout = window_days*m+h

# tt = length(algn$y)
nyearsin = 5 # number of years in sample
tin = 1:(nyearsin*12*m) #5 years insample
tout = (nyearsin*12*m)+h



##########################################################################################################################################################

###########################################################
#-------------- Pooled RUMIDAS model ---------------------#
###########################################################

######################################
#--------- Data preparation ---------#
######################################
#### Prepare data for (mixed-frequency) lags
## IF YOU USE THE POOLED RUMIDAS USE flexLF = T, FOR THE ORIGINAL RUMIDAS flexLF = F
# For pooled RUMIDAS
algn = align_data(x_hf = x_hf, y_lf = y_lf, p_hf = p_hf, p_lf = p_lf, m = m, h = h, flexLF = T)

# other input parameters
inp_mflags = inputs_algorithm(algn$k_hf, algn$k_lf, p_hf, p_lf, m, penalty = "HIER", flexLF = T)
group_index_mflags = inp_mflags$group_index
penalty_mflags = inp_mflags$penalty
columns_LF_mflags = inp_mflags$columns_LF
columns_HF_mflags = inp_mflags$columns_HF
# the inp_mflags$penalty_values returns are the recency-based penalty values: 
# however it is also possble to hard-code the penalty values if one has a specific penalty structure in mind
penalty_values_mflags = inp_mflags$penalty_values


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
# mflags-HIER: HIER on LF and HF
ymflags_hier_hats = Xmflags_outs %*% betas_mflags_hier[index_forecast_mflags]
ymflags_hier_hat = sigma_mflags_y * ymflags_hier_hats + mu_mflags_y

# error
ymflags_out - ymflags_hier_hat

######## TUNING PARAMETER SELECTION BY TIME SERIES CROSS-VALIDATION ##########
# Warning takes long 
# mflags_hier: hierarchical on LF and HF
mflags_hier_out = sparseRUMIDAS_tscv(y = ymflags_ins, X = Xmflags_ins, cvcut = 0.2*nrow(Xmflags_ins), 
                                     group_index = group_index_mflags, penalty = penalty_mflags, penalty_values = penalty_values_mflags, 
                                     columns_LF = columns_LF_mflags, columns_HF = columns_HF_mflags,
                                     l_lambda = l_lambda, post_lasso = TRUE, epsilon = 10^-4, max_iter = 500)
betas_mflags_hier = mflags_hier_out$betas
lambda_tscv_min = mflags_hier_out$lambda_min
CVscore = mflags_hier_out$MAFE_avg

##########################################################################################################################################################

############################################################
#-------------- Standard RUMIDAS model ---------------------#
############################################################
## IF YOU USE THE POOLED RUMIDAS USE flexLF = T, FOR THE ORIGINAL RUMIDAS flexLF = F
algn = align_data(x_hf = x_hf, y_lf = y_lf, p_hf = p_hf, p_lf = p_lf, m = m, h = h)
inp_mflags = inputs_algorithm(algn$k_hf, algn$k_lf, p_hf, p_lf, m, penalty = "HIER")
group_index_mflags = inp_mflags$group_index
penalty_mflags = inp_mflags$penalty
columns_LF_mflags = inp_mflags$columns_LF
columns_HF_mflags = inp_mflags$columns_HF
# the inp_mflags$penalty_values returns are the recency-based penalty values: 
penalty_values_mflags = inp_mflags$penalty_values

period_forecast = algn$hf_period_index[tout] 
datamat = prep_eqei(y = algn$y[tin], X = algn$X[tin,], hf_period_index = algn$hf_period_index[tin], mi = period_forecast, intercept = FALSE)

ymflags_in = datamat[,1]
Xmflags_in = datamat[,-1]
ymflags_out = algn$y[tout]
Xmflags_out = algn$X[tout,]

##### before penalization still standardize!! #####
Xmflags_ins = scale(Xmflags_in)
ymflags_ins = scale(ymflags_in)

mu_mflags_y = attributes(ymflags_ins)$`scaled:center`
sigma_mflags_y = attributes(ymflags_ins)$`scaled:scale`
mu_mflags_X = attributes(Xmflags_ins)$`scaled:center`
sigma_mflags_X = attributes(Xmflags_ins)$`scaled:scale`


############# with LF #######################
# mflags_hier: hierarchical on LF and HF
mflags_hier_out = sparseRUMIDAS(y = ymflags_ins, X = Xmflags_ins, group_index = group_index_mflags, penalty = penalty_mflags, penalty_values = penalty_values_mflags, 
                                columns_LF = columns_LF_mflags, columns_HF = columns_HF_mflags,
                                l_lambda = l_lambda, post_lasso = TRUE, epsilon = 10^-4, max_iter = 500)

betas_mflags_hier = mflags_hier_out$betas[,mflags_hier_out$bic_index]
lambda_bic_periodhf = mflags_hier_out$lambdaSeq[mflags_hier_out$bic_index]

# standardize X for HIER out-of-sample
Xmflags_outs = ( Xmflags_out - mu_mflags_X ) / sigma_mflags_X
ymflags_hier_hats = Xmflags_outs %*% betas_mflags_hier
ymflags_hier_hat = sigma_mflags_y * ymflags_hier_hats + mu_mflags_y

# error
ymflags_out - ymflags_hier_hat
