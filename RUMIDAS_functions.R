######## HAR model ########
prep_har <- function(RV, order = c(1,5,22), h = 1, leadlag = 0){
  data_lags = embed(RV, dimension = h+1)
  MA_1 = rollmean(data_lags[,ncol(data_lags)], k = order[2])
  MA_2 = rollmean(data_lags[,ncol(data_lags)], k = order[3])
  ydata = data_lags[-c(1:(nrow(data_lags)-length(MA_2))),1]
  Xdata = cbind(1,data_lags[-c(1:(nrow(data_lags)-length(MA_2))),ncol(data_lags)], MA_1[-c(1:(length(MA_1)-length(MA_2)))], MA_2)
  colnames(Xdata) = c(paste0("beta", c(0,order)))
  if(leadlag == 1){
    ydata = ydata[-c(1:max(order))]
    Xdata = Xdata[-c(1:max(order)),]
  }

  return(list("ydata" = ydata, "Xdata" = Xdata))
}

align_data <- function(x_hf, y_lf, p_hf, p_lf, m, h = 1, flexLF = FALSE, leadlag = 0){ 
  #!!! first variable of x_hf needs to be dependent variable
  # k_hf, k_lf: number of hf and lf variables
  # p_hf, p_lf: number of hf and lf lags
  # m: frequency mismatch 
  # h: forecast horizon
  # flexLF: if using the pooled RU-MIDAS set this to TRUE, for standard RUMIDAS = FALSE
  
  # leadlag = -1 # lead 
  # leadlag = 0 # usual p_lf
  # leadlag = 1 # lag more
  if(is.null(ncol(x_hf))){
    k_hf = 1
  }else{
    k_hf = ncol(x_hf)
  }
  if(is.null(ncol(y_lf))){
    k_lf = 1
  }else{
    k_lf = ncol(y_lf)
  }
  
  p = max(p_lf+leadlag, ceiling(p_hf/m))
  
  hf_lags = embed(x_hf, p*m+1)
  lf_lags = embed(y_lf, p+1)
  LF = matrix(rep(lf_lags[,c((k_lf*(leadlag+1)+1):(k_lf*(leadlag+1)+p_lf*k_lf))], each = m), nrow(lf_lags)*m , p_lf*k_lf)
  
  if(flexLF){
    LF_dummy = matrix(0, nrow(LF), ncol(LF)*m)
    HF_index = matrix(1:nrow(hf_lags), nrow(lf_lags), m, byrow = T)
    for(i in 1:m){
      LF_dummy[HF_index[,i],((k_lf*p_lf)*(i-1)+1):((k_lf*p_lf)*i)] = LF[HF_index[,i],]
    }
    LF = LF_dummy
  }
  
  HF = hf_lags[,c((k_hf+1):(k_hf+p_hf*k_hf))]
  X = cbind(LF, HF) # Lagged predictor matrix
  y = hf_lags[,1, drop = F]
  hf_period_index = rep(1:m, times = length(y)/m)
  
  if(h > 1){ #h-step forecast
    y = y[-(1:(h-1)), ,drop = F]
    X = X[-c((nrow(X)-h+2):nrow(X)),]
    hf_period_index = hf_period_index[-(1:(h-1))]
  }
  
  return(list("y"= y, "X" = X, "hf_period_index" = hf_period_index, "k_lf" = k_lf, "k_hf" = k_hf, "h" = h))
}

align_mfhar <- function(RV, y_lf, p_lf, order = c(1,5,22), h = 1, m = 22, flexLF  = F, leadlag = 0){
  # leadlag = -1 # lead 
  # leadlag = 0 # usual 
  # leadlag = 1 # lag more
  k_hf = 1
  x_hf = RV
  
  if(is.null(ncol(y_lf))){
    k_lf = 1
  }else{
    k_lf = ncol(y_lf)
  }
  
  if(p_lf > 22){
    stop("p_lf needs to be smaller (or equal) than 22")
  }
 
  hf_lags = embed(x_hf, 2)
  lf_lags = embed(y_lf, leadlag+p_lf+1)
  
  MA_1 = rollmean(hf_lags[,2], k = order[2])
  MA_2 = rollmean(hf_lags[,2], k = order[3])
  HF = cbind(hf_lags[-c(1:(nrow(hf_lags)-length(MA_2))),2], MA_1[-c(1:(length(MA_1)-length(MA_2)))], MA_2)
  y = hf_lags[-c(1:(nrow(hf_lags)-length(MA_2))),1, drop = F]
  
  LF = matrix(rep(lf_lags[,c((k_lf*(leadlag+1)+1):(k_lf*(leadlag+1)+p_lf*k_lf))], each = m), nrow(lf_lags)*m , p_lf*k_lf) #c((k_lf+1):(k_lf+p_lf*k_lf))
  
  if(flexLF){
    LF_dummy = matrix(0, nrow(LF), ncol(LF)*m)
    HF_index = matrix(1:nrow(LF), nrow(lf_lags), m, byrow = T)
    for(i in 1:m){
      LF_dummy[HF_index[,i],((k_lf*p_lf)*(i-1)+1):((k_lf*p_lf)*i)] = LF[HF_index[,i],]
    }
    LF = LF_dummy
  }
  
  if(leadlag == 1){
    Lag = 1
  }else{
    Lag = 0
  }
  if(leadlag == -1){
    LF = LF[-c(1:m),]
  }else{
    Lead = 0
  }
  
  
  ydata = y[((p_lf+Lag-1)*m+1):nrow(HF),,drop =F]
  Xdata = cbind(LF, HF[((p_lf+Lag-1)*m+1):nrow(HF),])
  colnames(Xdata) <- NULL 
  #colnames(Xdata) = c( paste0(c(paste0("LF", 1:k_lf)), rep(paste0("_t-",1:p_lf),each = k_lf)), paste0("beta", order))
  hf_period_index = rep(1:m, times = length(ydata)/m)
  
  if(h > 1){ #h-step forecast
    ydata = ydata[-(1:(h-1)), ,drop = F]
    Xdata = Xdata[-c((nrow(Xdata)-h+2):nrow(Xdata)),]
    hf_period_index = hf_period_index[-(1:(h-1))]
  }
  
  return(list("y" = ydata, "X" = Xdata, "hf_period_index" = hf_period_index, "k_lf" = k_lf, "h" = h))
}

# help function for standard RUMIDAS to do equation-by-equation estimation 
prep_eqei <- function(y, X, hf_period_index, mi, intercept = FALSE){
  
  n = nrow(X) #(nrow(quarterly) - p_lf)*m
  #n_lf = n/m  #nrow(quarterly) - p_lf
  mi_index = which(hf_period_index == mi)
  #HF_index = matrix(1:n, n_lf, m, byrow = T)
  
  if(intercept == FALSE){
    yi = y[mi_index]
    Xi = X[mi_index,]
    data_mat = cbind(yi, Xi)
  }
  if(intercept == TRUE){
    yi = y[mi_index]
    Xi = X[mi_index,]
    data_mat = cbind(yi, 1, Xi)
  }
  return(data_mat)
}


