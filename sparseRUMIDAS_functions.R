sparseRUMIDAS <- function(y, X, group_index, penalty, penalty_values,
                          columns_LF, columns_HF, 
                          l_lambda = 10, lambdas = NULL, lambda_weight = rep(1, max(group_index)), 
                          post_lasso = TRUE, epsilon = 10^-4, max_iter = 500){
  
  # Input: 
  # y: dependent high-frequency (hf) variable 
  # X: regressor matrix containing lags of low-frequency (lf) and (hf) variable
  # group_index: indicates the group structure (of columns in X)
  # penalty: type of penalty "HIER" (hierarchical) or "L1" (lasso)
  # penalty_values: penalty priority values that indicate the hierarchical structure (needed for hierarchical)
  # columns_LF, columns_HF: indicate the columns that belong to lf and hf variables in X
  # l_lambda: number of tuning parameters if no tuning grid is provided
  # lambdas: provide tuning parameter
  # lambda_weight: to vary tuning parameter for each regressor (optional)
  # post_lasso: whether the post lasso is applied (TRUE) or not (FALSE)
  # epsilon: a small positive numeric value giving the tolerance for convergence in the proximal gradient algorithm.
  # max_iter: scaler, specifies maximum number of iterations in proximal gradient algorithm
  
  # Output:
  # y: dependent high-frequency (hf) variable 
  # X: regressor matrix containing lags of low-frequency (lf) and (hf) variable
  # betas: Matrix of estimated coefficients across different lambdas
  #        coefficients in column 1 correspond to the sparsest solution, in the last column to the most dense solution
  # post_lasso: whether the post lasso is applied (TRUE) or not (FALSE)
  # columns_LF, columns_HF: indicate the columns that belong to lf and hf variables in X
  # lambdaSeq: tuning parameter grid
  # bic: value BIC (tuning parameter selection)
  # bic_index: index lambda with min BIC 
  # penalty: type of penalty "HIER" (hierarchical) or "L1" (lasso)
  # penalty_values: penalty priority values that indicate the hierarchical structure (needed for hierarchical)
  # iter: vector containing iterations until convergence for each lambda 

  
  ydata = matrix(y, ncol = 1)
  
  # Case 2: Lags of LF variables are penalized 
  Xdata = X
  Xfull = X # for output
  RUMIDASmodel = model_prox(y = ydata, X = Xdata)
  s = step_size(Z = Xdata)
  
  if(is.null(lambdas)){
    lambdas = lambda_line(y = ydata, X = Xdata, group_index = group_index, penalty = penalty, penalty_values = penalty_values,
                          l_lambda = l_lambda, lambda_weight = lambda_weight, multipl_mat = RUMIDASmodel, s = s, epsilon = epsilon, max_iter = max_iter)
    
    PLS_results <- PLS(y = ydata, X = Xdata, group_index = group_index, penalty = penalty, penalty_values = penalty_values,
                       lambda_seq = lambdas, lambda_weight = lambda_weight, epsilon = epsilon, max_iter = max_iter)
    
    betas_PLS = PLS_results$betas
    iter_PLS = PLS_results$iters
    
    if(post_lasso == TRUE){
      post = apply(betas_PLS, 2, FUN = postlasso, ydata = ydata, Xdata = Xdata)
      betas_PLS = post
    }
    
    # BIC selection
    bic = apply(betas_PLS, 2, BIC, ydata = ydata, Xdata = Xdata)
    bic_index = which.min(bic)
    
  }else{ #just one lambda
    PLS_results <- PLS(y = ydata, X = Xdata, group_index = group_index, penalty = penalty, penalty_values = penalty_values,
                       lambda_seq = lambdas, lambda_weight = lambda_weight, epsilon = epsilon, max_iter = max_iter)
    betas_PLS = PLS_results$betas
    iter_PLS = PLS_results$iters
    
    if(post_lasso == TRUE){
      post = postlasso(ydata, Xdata, betas_PLS)
      betas_PLS = post
    }
    bic = bic_index = NULL
  }
  
  out <- list("y" = ydata, "X" = Xfull,
              "betas" = betas_PLS, "post_lasso" = post_lasso,
              "columns_LF" = columns_LF, "columns_HF" = columns_HF,
              "lambdaSeq"= lambdas, "bic" = bic, "bic_index" = bic_index, 
              "penalty" = penalty,"penalty_values" = penalty_values,
              "iter"= iter_PLS)
  return(out) 
}


# prox for an entire tuning parameter sequence (using warm-starts)
PLS <- function(y, X, group_index, penalty, penalty_values, lambda_seq, lambda_weight = rep(1, max(group_index)),
                       epsilon = 10^-4, max_iter = 500){
  
  # Input
  # y: vector Nx1
  # X: independent variables NxK (can be either lagged predictor matrix or seasonality matrix with fourier terms)
  # penalty: type of penalty for X c("HIER", "L1") if X lagged matrix and "GL" if seasonality matrix
  # penalty_values: penalty values for penalty X (for "L1" rep(1, ncol(Z)), for "HIER" hierarchical penalty structure, for "GL" group structure)
  # lambda_seq: tuning parameter sequence
  # epsilon: convergence criterion (default = 10^-4)
  # max_iter: maximum number of iterations allowed (default = 500)
  # s: step size for proximal gradient algorithm 
  # multipl_mat: contains several matrices to speed up calculation of gradient (saves recalculating them each time)) (default NULL)
  # beta_warm: warm start for estimates (default NULL)
  
  # Output: 
  # betas: matrix Kxlength(lambda_seq) of estimates for each element in lambda sequence 
  # iters: number of iterations needed until convergence for each element in lambda grid 
  
  # Construct necessary data matrices
  multipl_mat = model_prox(y = y, X = X)
  
  s = step_size(Z = X)
  
  beta_warm = rep(0, ncol(X))
  
  # save betas and iterations
  betas = matrix(0, nrow = ncol(X), ncol = length(lambda_seq))
  rownames(betas) <- colnames(X)
  iters = rep(0, length(lambda_seq))
  
  for(i in 1:length(lambda_seq)){
    list_prox <- prox(y = y, X = X, group_index = group_index, penalty = penalty, penalty_values = penalty_values, 
                                    lambda = lambda_seq[i], lambda_weight = lambda_weight, 
                                    epsilon = epsilon, max_iter = max_iter, 
                                    s = s, 
                                    multipl_mat = multipl_mat, beta_warm = beta_warm)
    betas[,i] = list_prox$beta
    iters[i] = list_prox$iter
    
    # warm start for next, smaller lambda (previous solution is new starting value)
    beta_warm = betas[,i]
  }
  results = list("betas" = betas, "iters" = iters)
  return(results)
}


rumidas_standardize <- function(y, Xdummy, m){
  n = nrow(Xdummy) #(nrow(quarterly) - p_lf)*m
  n_lf = n/m  #nrow(quarterly) - p_lf
  K = ncol(Xdummy)
  k = K/m
  y = matrix(y, ncol = 1)
  
  # equation-by-equation (m equations, one equation for each i = 1,..,m)
  Xd_s = matrix(0, nrow(Xdummy), ncol(Xdummy))
  yd_s = y
  mu_y = rep(NA,m)
  mu_X = rep(NA,ncol(Xdummy))
  sigma_y = rep(NA,m)
  sigma_X = rep(NA,ncol(Xdummy))
  
  for(i in 1:m){
    mat = Xdummy[,(k*(i-1)+1):(k*i)]
    index = which(rowSums(mat != 0.0) > 0)
    Xis = scale(mat[index,])
    yis = scale(y[index,], center = TRUE, scale = TRUE)
    Xd_s[index,(k*(i-1)+1):(k*i)] = Xis
    yd_s[index,] = yis
    
    mu_y[i] = attributes(yis)$`scaled:center`
    sigma_y[i] = attributes(yis)$`scaled:scale`
    mu_X[(k*(i-1)+1):(k*i)] = attributes(Xis)$`scaled:center`
    sigma_X[(k*(i-1)+1):(k*i)] = attributes(Xis)$`scaled:scale`
    
  }
  
  Xd_s = Matrix(Xd_s, sparse = TRUE) 
  
  out = list("yd_s" = yd_s, "Xd_s" = Xd_s, "mu_y" = mu_y, "sigma_y" = sigma_y, "mu_X" = mu_X, "sigma_X" = sigma_X)
  return(out)
}

rumidas_standardize_flexLF <- function(y, Xdummy_LF, X_HF, m){
  n = nrow(Xdummy_LF) #(nrow(quarterly) - p_lf)*m
  n_lf = n/m  #nrow(quarterly) - p_lf
  K = ncol(Xdummy_LF)
  k = K/m
  y = matrix(y, ncol = 1)
  
  ######### HF variable
  X_HF_ins = scale(X_HF)
  y_ins = scale(y)
  
  mu_y = attributes(y_ins)$`scaled:center`
  sigma_y = attributes(y_ins)$`scaled:scale`
  mu_HF_X = attributes(X_HF_ins)$`scaled:center`
  sigma_HF_X = attributes(X_HF_ins)$`scaled:scale`
  
  ####### LF Variable equation-by-equation (m equations, one equation for each i = 1,..,m)
  Xd_s = matrix(0, nrow(Xdummy_LF), ncol(Xdummy_LF))
  mu_X = rep(NA,ncol(Xdummy_LF))
  sigma_X = rep(NA,ncol(Xdummy_LF))
  
  for(i in 1:m){
    mat = Xdummy_LF[,(k*(i-1)+1):(k*i), drop = F]
    index = which(rowSums(mat != 0.0) > 0)
    Xis = scale(mat[index,])
    Xd_s[index,(k*(i-1)+1):(k*i)] = Xis

    mu_X[(k*(i-1)+1):(k*i)] = attributes(Xis)$`scaled:center`
    sigma_X[(k*(i-1)+1):(k*i)] = attributes(Xis)$`scaled:scale`
  }
  
  X_s_out = cbind(Xd_s, X_HF_ins)
  mu_X_out = c(mu_X, mu_HF_X)
  sigma_X_out = c(sigma_X, sigma_HF_X)
  
  #Xd_s = Matrix(Xd_s, sparse = TRUE) 
  
  out = list("yd_s" = y_ins, "Xd_s" = X_s_out, "mu_y" = mu_y, "sigma_y" = sigma_y, "mu_X" = mu_X_out, "sigma_X" = sigma_X_out)
  return(out)
}



# prepare data for algorithm
inputs_algorithm <- function(k_hf, k_lf, p_hf, p_lf, m, penalty = c("L1", "HIER"), flexLF = F){ #intercept = TRUE )
  #!!! first variable of x_hf needs to be dependent variable !!!
  # k_hf, k_lf: number of hf and lf variables
  # p_hf, p_lf: number of hf and lf lags
  # m: frequency mismatch 
  # penalty: L1 (lasso) or HIER (hierarchical)
  # flexLF: if using the pooled RU-MIDAS set this to TRUE, for standard RUMIDAS = FALSE
 
  # Group index
  group_index = groups_function(k_hf, p_hf, k_lf, p_lf, m, flexLF = flexLF)

  
  # Penalty values
  if(penalty == "L1"){
    if(flexLF){
      penalty_values = rep(1, (k_lf*p_lf*m + k_hf*p_hf))
    }else{
      penalty_values = rep(1, (k_lf*p_lf + k_hf*p_hf))
    }
  }
  if(penalty == "HIER"){
    if(flexLF){
      penalty_values_LF = rep(c(matrix(c(1:(m*p_lf)), ncol = m, byrow = T)), each = k_lf)
    }else{
      penalty_values_LF = rep(1:p_lf, each = k_lf)
    }
    penalty_values_HF = rep(1:p_hf, each = k_hf)
    penalty_values = c(penalty_values_LF, penalty_values_HF) #rep(c(penalty_values_LF, penalty_values_HF), m)
  }
  
  # Columns LF / HF
  if(k_lf == 0){
    columns_LF = NA
    columns_HF = c((k_lf*p_lf+1):(k_lf*p_lf+k_hf*p_hf))
  }else{
    if(flexLF){
      columns_LF = c(1:(k_lf*p_lf*m))
      columns_HF = c((k_lf*p_lf*m+1):(k_lf*p_lf*m+k_hf*p_hf))
    }else{
      columns_LF = c(1:(k_lf*p_lf))
      columns_HF = c((k_lf*p_lf+1):(k_lf*p_lf+k_hf*p_hf))
    }
  }
    
  out <- list("group_index" = group_index, "penalty" = penalty, "penalty_values" = penalty_values, "columns_LF" = columns_LF, "columns_HF" = columns_HF)
  
  return(out)
}


groups_function <- function(k_hf, p_hf, k_lf, p_lf, m, flexLF = FALSE){
  groups_LF = rep(1:k_lf, times = p_lf)
  if(flexLF){
    groups_LF = rep(groups_LF, m)
  }
  groups_HF = rep((k_lf+1):(k_lf+k_hf), times = p_hf)
  groups_m1 = c(groups_LF, groups_HF)
  
  group_index = groups_m1
  return(group_index)
}


# Accelerated proximal gradient method for one penalized group
prox <- function(y, X, group_index, penalty, penalty_values, lambda, lambda_weight = rep(1, max(group_index)), 
                 epsilon = 10^-4, max_iter = 500, 
                 s = NULL, 
                 multipl_mat = NULL, beta_warm = NULL){
  
  # Input
  # y: vector Nx1
  # X: independent variables NxK (can be either lagged predictor matrix or seasonality matrix with fourier terms)
  # penalty: type of penalty for X c("HIER", "L1") if X lagged matrix and "GL" if seasonality matrix
  # penalty_values: penalty values for penalty X (for "L1" rep(1, ncol(Z)), for "HIER" hierarchical penalty structure, for "GL" group structure)
  # lambda: tuning parameter
  # epsilon: convergence criterion (default = 10^-4)
  # max_iter: maximum number of iterations allowed (default = 500)
  # s: step size for proximal gradient algorithm 
  # multipl_mat: contains several matrices to speed up calculation of gradient (saves recalculating them each time)) (default NULL)
  # beta_warm: warm start for estimates (default NULL)
  
  # Output: 
  # beta: estimates 
  # iter: number of iterations needed until convergence
  
  
  # Construct necessary data matrices for gradient
  if(is.null(multipl_mat)){
    multipl_mat = model_prox(y = y, X = X)
  }
  y_X = multipl_mat$y_X
  XtranspX = multipl_mat$XtranspX
  
  # step size
  if(is.null(s)){
    s = step_size(Z = X)
  }
  
  # beta = 0 or beta_warm
  if(!is.null(beta_warm) & length(beta_warm) == ncol(X)){
    beta_oldold <- beta_old <- beta_warm[c(1:ncol(X))]
  }else{
    beta_oldold <- beta_old <- rep(0, ncol(X))
  }
  
  r <- 2 # iterations
  convergence <- FALSE
  
  while(convergence == FALSE & r < max_iter){
    r <- r+1
    beta_copy <- beta_old # make copy to compare convergence later 
    
    for(g in 1:max(group_index)){
      index = which(group_index == g)
      
      beta_hat <- beta_old[index] + ((r-2)/(r+1))*(beta_old[index]-beta_oldold[index])
      beta_oldold[index] <- beta_old[index]
      beta_old[index] <- beta_hat
      
      gradient_beta = -(y_X - t(beta_old) %*% XtranspX)[,index] 
      #gradient_beta = -(y_X[index] - t(beta_old) %*% XtranspX[,index])  #CHECK THIS!!

      update_beta = beta_hat - s*gradient_beta
      
      if(penalty == "L1"){
        # Update elementwise soft thresholding 
        beta_new <- softelem(c(update_beta), lambda = s*lambda*lambda_weight[g])  
      }
      if(penalty == "HIER"){
        # Update groupwise soft thresholding (hierarchical)
        beta_new <- softnestedgroup(c(update_beta), lambda = s*lambda*lambda_weight[g], group_indices = penalty_values[index])  
      }
      # Save new updated values 
      beta_old[index] <- beta_new
    }
    # Check if estimator has reached convergence
    if(max(abs(beta_copy-beta_old))<=epsilon){
      convergence = TRUE
    }
  }
  return(list("beta" = beta_old, "iter" = r))
 }  


# data matrix inputs for prox
model_prox <- function(y, X){
  # Input
  # y: vector Nx1
  # X: transformed X matrix, NxP matrix that contains eigenvectors of group if penalty "HIER" else Z = X
  
  # Output: 
  # y_X: 1xP matrix (useful for calculation of gradient (saves recalculating them each time))
  # XtranspX: PxP matrix (useful for calculation of gradient (saves recalculating them each time))
  # s: step size for proximal gradient algorithm 
  
  # Save values for calculation of gradient (saves recalculating them each time)
  y_X = as.matrix(t(y)%*%X)
  XtranspX = crossprod(X)
  
  out<-list("y" = y, "X" = X, "y_X" = y_X, "XtranspX" = XtranspX)
  return(out)
}


# step size for proximal gradient algorithm
step_size <- function(Z,X = NULL){
  ZX = cbind(Z,X)
  if(ncol(ZX)>2){
    ZXsvd = svds(ZX, k = 1, nu = 0, nv = 0)$d
  }else{
    ZXsvd = svd(ZX, nu = 0, nv = 0)$d[1]
  }
  s = (ZXsvd)^(-2)
  
  return(s)
}

# elementwise soft thresholding
softelem <- function(x, lambda){
  r <- sign(x)*pmax(abs(x)-lambda,0)
  return(r)
}


# groupwise soft thresholding
softgroup <-function(x, lambda, group_indices){
  r <- rep(NA, length(x))
  
  for(h in 1:max(group_indices)){
    subset <- which(group_indices==h)
    r_gh <- x[subset]
    r_gh <- max(0,1-lambda/norm(r_gh,"2"))*r_gh
    r[subset] <- r_gh
  }
  return(r)
}


# groupwise soft thresholding with nested groups
softnestedgroup <- function(x, lambda, group_indices){
  r <- x
  total_length = length(r)
  for(h in 1:max(group_indices)){
    if(h==1){
      r_gh = r
      r_gh <- max(0,1-lambda/norm(r_gh,"2"))*r_gh
      r = r_gh
      indices_deleted_groups <- c()
    }
    if(h!=1){
      indices_deleted_groups <- c(indices_deleted_groups, which(group_indices==(h-1)))
      
      # Take only elements of variable r that are still allowed to updated
      r_gh <- r[-indices_deleted_groups]
      # weight w
      #group_length = length(r_gh)
      #w = ((total_length - group_length) + 1) 
      # Do groupwise soft thresholding
      r_gh <- max(0,1-(lambda)/norm(r_gh,"2"))*r_gh 
      # Update the values in r (and keep unchanged for future iterations) which belong to group h 
      r[-indices_deleted_groups] = r_gh
    }
  }
  return(r)
}

# find lambda max for prox
lambdamax_search <- function(y, X, group_index, penalty, penalty_values, lstart, lambda_weight = rep(1, max(group_index)), 
                              multipl_mat = NULL, s = NULL, 
                              epsilon = 10^-4, max_iter = 500){
  
  if(is.null(multipl_mat)){
    multipl_mat = model_prox(y = y, X = X)
  }
  if(is.null(s)){
    s = step_size(Z = X)
  }
  
  lambda_left = 0
  lambda_right = lstart
  thresh = 10 
  while(thresh > 0.0001){
    result <- prox(y = y, X = X, group_index = group_index, penalty = penalty, penalty_values = penalty_values, 
                          lambda = lambda_right, lambda_weight = lambda_weight, 
                          epsilon = epsilon, max_iter = max_iter, 
                          s = s, multipl_mat = multipl_mat)
    
    if(min(lambda_weight == 1)){ # all betas are penalized
      if(max(abs(result$beta)) == 0){
        lambda_works = lambda_right
        lambda_right = (lambda_left+lambda_right)/2
      }else{
        lambda_left = lambda_right
        lambda_right = lambda_right*1.5
      }
    }else{
      g = which(lambda_weight!=0) #subset of betas are penalized
      index = which(group_index %in% g)
      if(max(abs(result$beta[index])) == 0){
        lambda_works = lambda_right
        lambda_right = (lambda_left+lambda_right)/2
      }else{
        lambda_left = lambda_right
        lambda_right = lambda_right*1.5
      }
    }

    thresh = abs(lambda_right - lambda_left)
  }
  return(lambda_works)
}


# Construction of tuning parameter sequence (1-dimensional)
lambda_line <- function(y, X, group_index, penalty, penalty_values, l_lambda, lambda_weight = rep(1, max(group_index)),
                        multipl_mat = NULL, s = NULL, 
                        epsilon = 10^-4, max_iter = 500){
  
  if(is.null(multipl_mat)){
    multipl_mat = model_prox(y = y, X = X)
  }
  if(is.null(s)){
    s = step_size(Z = X)
  }
  
  matX = t(X)%*%y
  lambda_max_start_X = max(abs(matX))
  
  lambda_max = lambdamax_search(y = y, X = X, group_index = group_index, penalty = penalty, penalty_values = penalty_values,
                                 lstart = lambda_max_start_X, lambda_weight = lambda_weight, 
                                 multipl_mat = multipl_mat, s = s,
                                 epsilon = epsilon, max_iter = max_iter)
  
  lambda_min = lambda_max/10^4
  lambda_seq <- c(exp(seq(log(lambda_max),log(lambda_min), length = l_lambda)))
  
  return(lambda_seq)
}

# Info criteria to select tuning parameter
BIC <- function(ydata, Xdata, beta){
  # Input:
  # ydata: y N x 1 vector
  # Xdata: n x P matrix
  # beta: P vector
  # ic: type of information criterion ("AIC" or "BIC")
  
  # Output:
  # Information criterion ("AIC" or "BIC")
  if(is.null(dim(Xdata))){
    n = length(Xdata)
  }else{
    n = nrow(Xdata)
  }
  yhat <- Xdata%*%beta
  fit <- sum((ydata-yhat)^2)/n 
  dF <- sum(beta!=0) #degrees of freedom = number of non-zero coefficient
  
  bic <- log(fit) + (1/n)*dF*log(n)
  return(bic)
}

postlasso <- function(ydata, Xdata, betastack){
  beta_post <- rep(0, length(betastack))
  
  subset <- which(betastack!=0)
    
  if(length(subset)!=0){ # at least one non-zero
    Xsub <- Xdata[,subset]
    OLS_post <- try(solve(crossprod(Xsub))%*%t(Xsub)%*%ydata, silent = TRUE)
    if(class(OLS_post)[1]!="try-error"){
      beta_post[subset] <- OLS_post
    }else{
      beta_post <- rep(NA, length(betastack))
    }
  } 
  return(beta_post)
}

sparseRUMIDAS_tscv <- function(y, X, cvcut = nrow(X)*0.2, group_index, penalty, penalty_values,
                               columns_LF, columns_HF, 
                               l_lambda = 10, lambda_choice = "min", lambda_weight = rep(1, max(group_index)), 
                               post_lasso = TRUE, epsilon = 10^-4, max_iter = 500){
  # Input: 
  # y: dependent high-frequency (hf) variable 
  # X: regressor matrix containing lags of low-frequency (lf) and (hf) variable
  # cvcut: number of observations used for forecast evaluation in the time series cross-validation procedure. The remainder is used for model estimation.
  # group_index: indicates the group structure (of columns in X)
  # penalty: type of penalty "HIER" (hierarchical) or "L1" (lasso)
  # penalty_values: penalty priority values that indicate the hierarchical structure (needed for hierarchical)
  # columns_LF, columns_HF: indicate the columns that belong to lf and hf variables in X
  # l_lambda: number of tuning parameters if no tuning grid is provided
  # lambda_choice = "min" or "1se"
  # lambda_weight: to vary tuning parameter for each regressor (optional)
  # post_lasso: whether the post lasso is applied (TRUE) or not (FALSE)
  # epsilon: a small positive numeric value giving the tolerance for convergence in the proximal gradient algorithm.
  # max_iter: scaler, specifies maximum number of iterations in proximal gradient algorithm
  
  # Output:
  # y: dependent high-frequency (hf) variable 
  # X: regressor matrix containing lags of low-frequency (lf) and (hf) variable
  # betas: Matrix of estimated coefficients across different lambdas
  #        coefficients in column 1 correspond to the sparsest solution, in the last column to the most dense solution
  # post_lasso: whether the post lasso is applied (TRUE) or not (FALSE)
  # columns_LF, columns_HF: indicate the columns that belong to lf and hf variables in X
  # lambdaSeq: tuning parameter grid
  # lambda_min: lambda with minimum MAFE (tuning parameter selection)
  # lambda_1se: lambda with minimum MAFE+1se (tuning parameter selection)
  # penalty: type of penalty "HIER" (hierarchical) or "L1" (lasso)
  # penalty_values: penalty priority values that indicate the hierarchical structure (needed for hierarchical)
  # iter: vector containing iterations until convergence for each lambda 
  # MAFE_avg: CV score (average) for each lambda
  # MAFE_matrix: CV scopre for each lambda for each iteration. 

  t <- nrow(X)
  n1 <- t-cvcut
  
  yfull = matrix(y, ncol = 1)
  ydata = matrix(y, ncol = 1)
  
  MAFE_matrix = matrix(NA, nrow = cvcut, ncol = l_lambda)
  rownames(MAFE_matrix) = paste0("n", 1:cvcut)
  colnames(MAFE_matrix) = paste0("lambdas", 1:l_lambda)
  
  # Case 2: Lags of LF variables are penalized (not going to programm the other case)
  
  Xdata = X
  Xfull = X # for output
  RUMIDASmodel = model_prox(y = ydata, X = Xdata)
  s = step_size(Z = Xdata)
  
  # Fit once on insample + validation set to find lambda sequence 
  lambdas = lambda_line(y = ydata, X = Xdata, group_index = group_index, penalty = penalty, penalty_values = penalty_values,
                        l_lambda = l_lambda, lambda_weight = lambda_weight, multipl_mat = RUMIDASmodel, s = s, epsilon = epsilon, max_iter = max_iter)
  
  for(n in 1:cvcut){
    ytrain = ydata[(n:(n1+n-1)),, drop = F]
    Xtrain = Xdata[(n:(n1+n-1)),]
    
    ytest = ydata[n1+n, drop = F]
    Xtest = Xdata[n1+n,]
    
    PLS_results_cv <- PLS(y = ytrain, X = Xtrain, group_index = group_index, penalty = penalty, penalty_values = penalty_values,
                          lambda_seq = lambdas, lambda_weight = lambda_weight, epsilon = epsilon, max_iter = max_iter)
    
    betas_PLS_cv = PLS_results_cv$betas
    
    if(post_lasso == TRUE){
      post = apply(betas_PLS_cv, 2, FUN = postlasso, ydata = ydata, Xdata = Xdata)
      betas_PLS_cv = post
    }
    
    # CV score
    ae_errors <- get_errors(beta_mat = betas_PLS_cv, ydata = ytest, Xdata = Xtest) # dimension 1 x l_lambda
    MAFE_matrix[n,] = ae_errors
  }
  
  CV_score = colMeans(MAFE_matrix)
  min_CV = which.min(CV_score)
  se_min_CV = sd(MAFE_matrix[,min_CV])/sqrt(cvcut)
  SE_rule = CV_score[min_CV]+se_min_CV
  sparse_CV = which(CV_score <= SE_rule, arr.ind = TRUE)[1]
  
  lambda_min <- lambdas[min_CV]
  lambda_1se <- lambdas[sparse_CV]
  
  
  if(lambda_choice == "min"){
    PLS_results <- PLS(y = ydata, X = Xdata, group_index = group_index, penalty = penalty, penalty_values = penalty_values,
                       lambda_seq = lambda_min, lambda_weight = lambda_weight, epsilon = epsilon, max_iter = max_iter)
  }
  if(lambda_choice == "1se"){
    PLS_results <- PLS(y = ydata, X = Xdata, group_index = group_index, penalty = penalty, penalty_values = penalty_values,
                       lambda_seq = lambda_1se, lambda_weight = lambda_weight, epsilon = epsilon, max_iter = max_iter)
  }
  
  betas_PLS = PLS_results$betas
  iter_PLS = PLS_results$iters
  
  if(post_lasso == TRUE){
    post = apply(betas_PLS, 2, FUN = postlasso, ydata = ydata, Xdata = Xdata)
    betas_PLS = post
  }
  
  out <- list("y" = ydata, "X" = Xfull,
              "betas" = betas_PLS, 
              "post_lasso" = post_lasso,
              "columns_LF" = columns_LF, "columns_HF" = columns_HF,
              "lambdaSeq"= lambdas, "lambda_min"=lambda_min,"lambda_1se"=lambda_1se, "lambda_min_index"=min_CV,"lambda_1se_index"=sparse_CV, 
              "penalty" = penalty,"penalty_values" = penalty_values,
              "iter"= iter_PLS, "MAFE_avg" = CV_score, "MAFE_all" = MAFE_matrix)
}


# Fitted values 
yhats_function <- function(beta_vec, Xdata){
  yhat <- Xdata%*%beta_vec
  return(yhat) # returns vectorized version
}

# Errors for multiple columns of betas
get_errors <- function(beta_mat, ydata, Xdata){
  yhats <- apply(beta_mat, 2, yhats_function, Xdata = Xdata)
  if(is.null(dim(yhats))){
    yhats = matrix(yhats, nrow = 1)
  }
  errors = apply(-yhats, 2, FUN = "+", ydata)
  if(is.null(dim(errors))){
    errors = matrix(errors, nrow = 1)
  }
  se <- apply(errors, 2, FUN = function(x){x^2})
  ae <- apply(errors, 2, FUN = function(x){abs(x)})
  return(ae) #absolute errors in paper
}

