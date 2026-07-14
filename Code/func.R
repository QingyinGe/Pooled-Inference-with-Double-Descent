getSepTestL2 = function(df_train_sep, df_test_sep, method){
  mean_train = colMeans(df_train_sep)
  
  if(method == "S:Linear"){
    cov_est = linshrink_cov(df_train_sep)
  }
  
  if(method == "S:Cor"){
    cov_est = getConstCorCov(df_train_sep)
  }
  
  if(method == "S:Var"){
    cov_est = diag(diag(cov(df_train_sep)) + 1e-10) #some experts have 0 variance
  }
  
  w = rowSums(solve(cov_est))
  w = w/sum(w)
  
  df_test_sep_demean = t(t(df_test_sep) - mean_train)
  
  err_vec = df_test_sep_demean %*% w
  
  mean(err_vec^2)
}

prepareFluData = function(df, which_season){
  df_1season = df %>% filter(Season == which_season) %>% dplyr::select(location, model_name, Model.Week, err) %>% 
    arrange(Model.Week) %>% pivot_wider(names_from = Model.Week, values_from = err) %>% arrange(location, model_name)
  
  colnames(df_1season)[3:ncol(df_1season)] = paste0(which_season, "/", colnames(df_1season)[3:ncol(df_1season)])
  
  df_1season
}

getRidge0TestErr = function(X_train, Y_train, X_test, Y_test, stdize = TRUE, n_experts = 17){
  n_train = nrow(X_train)
  n_prods = ncol(X_train)/(n_experts - 1)
  
  X_mu_train = colMeans(X_train)
  y_mu_train = colMeans(Y_train)
  
  if(stdize == TRUE){
    X_sd_train = apply(X_train, 2, sd)
    X_sd_train = replaceZeroWOne(X_sd_train) #some standard deviations are exactly 0
  }else if(stdize == "ProdConst"){
    X_sd_train = apply(X_train, 2, sd)
    X_sd_train = replaceZeroWOne(X_sd_train)
    benchmark_sd = rep(X_sd_train[1:(n_experts-1)], n_prods)
    ratio_train = matrix(X_sd_train/benchmark_sd, nrow = n_experts -1)
    X_sd_train = rep(colMeans(ratio_train), each = n_experts - 1)
  }else if(stdize == "ProdVar"){
    X_sd_train = apply(X_train, 2, sd)
    X_sd_train = replaceZeroWOne(X_sd_train)
    benchmark_sd = rep(X_sd_train[1:(n_experts-1)], n_prods)
    X_sd_train = X_sd_train/benchmark_sd
  }else{
    X_sd_train = rep(1, ncol(X_train))
  }
  
  X_train_scaled = t((t(X_train)-X_mu_train)/X_sd_train)
  Y_train_demean = t((t(Y_train)-y_mu_train))
  
  #Try to replace the following with X^T(XX^T)^{-1}Y? The demean part is non-trivial to handle.
  
  svd_train = svd(X_train_scaled)
  non_zero_idx = 1:(n_train-1)
  
  VDinv = svd_train$v[,non_zero_idx] %*% diag(1/svd_train$d[non_zero_idx])
  UTY = t(svd_train$u[,non_zero_idx]) %*% Y_train_demean #demean or without demean generates the same solution
  
  Y_pred = t((t(X_test) - X_mu_train)/X_sd_train) %*% VDinv %*% UTY
  
  err_test = Y_test - t(t(Y_pred) + y_mu_train)
  
  t(as.matrix(err_test))
}

getRidge0TestErrSim = function(X_train, Y_train, X_test, Y_test, stdize = TRUE, n_experts = 17){
  n_train = nrow(X_train)
  n_prods = ncol(X_train)/(n_experts - 1)
  
  X_mu_train = colMeans(X_train)
  y_mu_train = colMeans(Y_train)
  
  if(stdize == TRUE){
    X_sd_train = apply(X_train, 2, sd)
    X_sd_train = replaceZeroWOne(X_sd_train) #some standard deviations are exactly 0
  }else if(stdize == "ProdConst"){
    X_sd_train = apply(X_train, 2, sd)
    X_sd_train = replaceZeroWOne(X_sd_train)
    benchmark_sd = rep(X_sd_train[1:(n_experts-1)], n_prods)
    ratio_train = matrix(X_sd_train/benchmark_sd, nrow = n_experts -1)
    X_sd_train = rep(colMeans(ratio_train), each = n_experts - 1)
  }else if(stdize == "ProdVar"){
    X_sd_train = apply(X_train, 2, sd)
    X_sd_train = replaceZeroWOne(X_sd_train)
    benchmark_sd = rep(X_sd_train[1:(n_experts-1)], n_prods)
    X_sd_train = X_sd_train/benchmark_sd
  }else{
    X_sd_train = rep(1, ncol(X_train))
  }
  
  X_train_scaled = t((t(X_train)-X_mu_train)/X_sd_train)
  
  svd_train = svd(X_train_scaled)
  non_zero_idx = 1:(n_train-1)
  
  VDinv = svd_train$v[,non_zero_idx] %*% diag(1/svd_train$d[non_zero_idx])
  UTY = t(svd_train$u[,non_zero_idx]) %*% Y_train
  
  Y_pred = t((t(X_test) - X_mu_train)/X_sd_train) %*% VDinv %*% UTY
  
  err_test = Y_test - t(t(Y_pred) + y_mu_train)
  
  t(as.matrix(err_test))
}

getPCRTestErr = function(X_train, Y_train, X_test, Y_test, threshold = 0.8){
  n_train = nrow(Y_train)
  n_variable = ncol(Y_train)
  n_covariates = ncol(X_train)
  n_test = nrow(Y_test)
  
  if (n_train > n_covariates) {
    stop("it is not in the modern regime")
  }
  
  X_mu_train = colMeans(X_train)
  y_mu_train = colMeans(Y_train)
  
  X_sd_train = apply(X_train, 2, sd)
  X_sd_train = replaceZeroWOne(X_sd_train) #some standard deviations are exactly 0
  
  X_train_scaled = t((t(X_train)-X_mu_train)/X_sd_train)
  Y_train_demean = t((t(Y_train)-y_mu_train))
  
  svd_train = svd(X_train_scaled)
  pct_sum = cumsum(svd_train$d^2)/sum(svd_train$d^2)
  PCA_idx = 1:min(which(pct_sum > threshold))
  # non_zero_idx = 1:(n_train-1)
  
  if(length(PCA_idx) >1){
    VDinv = svd_train$v[,PCA_idx,drop = FALSE] %*% diag(1/svd_train$d[PCA_idx])
  }else{
    VDinv = svd_train$v[,PCA_idx,drop = FALSE]/svd_train$d[PCA_idx]
  }
  
  UTY = t(svd_train$u[,PCA_idx]) %*% Y_train_demean #demean or without demean generates the same solution
  
  # VDinv %*% UTY could be replaced by t(X_train_scaled) %*% ginv(X_train_scaled %*% t(X_train_scaled)) %*% Y_train_demean
  
  Y_pred = t((t(X_test) - X_mu_train)/X_sd_train) %*% VDinv %*% UTY
  
  err_test = Y_test - t(t(Y_pred) + y_mu_train)
  
  t(as.matrix(err_test))
}


replaceZeroWOne = function(x_vec){
  x_vec[which(x_vec == 0)] = 1
  x_vec
}


getGamma0 = function(n_experts){
  root_mat = diag(n_experts) - 1/n_experts #1/n_experts = 1/n_experts * one_vec %*% t(one_vec)
  
  Gamma0_raw = root_mat[, 1:(n_experts - 1)]
  A = t(Gamma0_raw) %*% Gamma0_raw
  eig = eigen(A)
  A_inv_half =
    eig$vectors %*%
    diag(1 / sqrt(eig$values)) %*%
    t(eig$vectors)
  
  Gamma0 = Gamma0_raw %*% A_inv_half
}

getTimeSeriesNames = function(data_freq, rank_idx){
  df = read.csv("CleanedData/M4ScaledData/M4-info.csv")
  
  df_sub = df %>% filter(SP == data_freq)
                             
  table_start_time = table(df_sub$StartingDate)
  
  chosen_start_time = names(sort(table_start_time, decreasing = TRUE))[rank_idx]
  
  df_chosen = df_sub %>% filter(StartingDate == chosen_start_time)
  
  df_chosen$M4id
}

getTimeSeriesCats = function(data_freq, rank_idx){
  df = read.csv("CleanedData/M4ScaledData/M4-info.csv")
  
  df_sub = df %>% filter(SP == data_freq)
  
  table_start_time = table(df_sub$StartingDate)
  
  chosen_start_time = names(sort(table_start_time, decreasing = TRUE))[rank_idx]
  
  df_chosen = df_sub %>% filter(StartingDate == chosen_start_time)
  
  round(table(df_chosen$category)/nrow(df_chosen), digits = 2)
}

sumEachProd = function(u_vec, n_products, n_experts){
  #speed up large matrix multiplication involving matrix E (Proposition 3)
  map_dbl(1:n_products,~sum(u_vec[1:n_experts + (.x-1) * n_experts ]))
}

getLambdaFromLinShrink = function(svd_X, X){ 
  # only use this function when n_products * n_experts >> n_obs
  # it avoids calculating the sample covariance matrix, which could be massive when n_products is large
  # We rescale linshrink_cov to be (n_obs-1) S + lambda I 
  
  n = length(svd_X$d)
  p = nrow(svd_X$v)
  
  m = sum(svd_X$d^2)/(n-1)/p
  
  d2 = (sum(svd_X$d^4)/(n-1)^2 - 2 * m * sum(svd_X$d^2)/(n-1) + m^2 *p)/p
  
  term1 = sum(apply(X, 1, function(x) sum(x^2)^2))
  
  term2 = sum(apply(X, 1, function(x) sum(svd_X$d[-n]^2 * (x %*% svd_X$v)^2)))/(n-1)
  
  term3 = n * sum(svd_X$d^4)/(n-1)^2
  
  b_bar2 = (term1 - 2* term2 + term3)/(n-1)^2/p
  
  b2 = min(d2, b_bar2)
  a2 = d2 - b2
  
  return(b2 * m * (n - 1)/a2)
}


getFastPoolLinearTestU = function(X, X_test, n_products, n_experts, n_obs){
  # We use the Woodbury inequality to derive this algorithm
  # We have checked that it returns the same result as using getPoolWeights
  # X_test has been debiased using train average (mu_vec or mean_in)
  X = scale(X, scale = FALSE) #demean
  svd_X = svd(X, nu = 0, nv = n_obs - 1)
  
  lin_lambda = getLambdaFromLinShrink(svd_X, X)
  
  EtV = apply(svd_X$v, 2, sumEachProd, n_products = n_products, n_experts = n_experts)
  EtU = apply(t(X_test), 2, sumEachProd, n_products = n_products, n_experts = n_experts)
  VtU = t(svd_X$v) %*% t(X_test)
  inv_EtinvE = 1/n_experts *
    (diag(n_products) + EtV %*%
       solve(n_experts * diag(lin_lambda/svd_X$d[-n_obs]^2 + 1) -  t(EtV) %*% EtV) %*% t(EtV))
  
  test = inv_EtinvE %*% (EtU - EtV %*% solve(diag(lin_lambda/svd_X$d[-n_obs]^2 + 1)) %*% VtU)
  
  return(test)
}


getPoolWeights = function(cov_inv, n_products, n_experts){
  # Formula from Proposition 3 
  E = diag(n_products) %x% matrix(1, n_experts, 1)
  weights = cov_inv %*% E %*% solve(t(E) %*% cov_inv %*% E)
  return(weights)
}



getPooledUFast = function(pool_idx, n_prods, n_experts, u_comb, true_data){
  #Leave-one-out result using map_dfc
  u_comb_agg = map_dfc(1:ncol(true_data), ~ 
                         getFastPoolLinearTestU(t(u_comb[pool_idx,-.x]), 
                                                t(as.matrix(u_comb[pool_idx,.x] - rowMeans(u_comb[pool_idx,-.x]))),
                                                n_products = n_prods, n_experts = n_experts, n_obs = ncol(true_data) - 1))
  
  u_comb_agg
}


getPooledU = function(idx, cov_type, m, n,f_comb, u_comb, true_data){
  #Leave-one-out result using map_dfc
  f_comb_agg = map_dfc(1:ncol(true_data), ~ aggForecast(f_comb[idx,.x], u_comb[idx,],
                                                        .x, m = m, n = n, n_days - 1,
                                                        cov_type = cov_type))
  
  u_comb_agg = f_comb_agg - true_data
  
  u_comb_agg
}

getSepUMethods = function(idx, methods, n, f_sep, u_sep, true_sep){
  n_methods = length(methods)
  u_comb_mat = matrix(NA, n_methods, ncol(u_sep))
  
  for(i in 1:n_methods){
    f_comb_agg = map_dfc(1:ncol(u_sep), ~ aggForecast(f_sep[idx,.x], u_sep[idx,],
                                                      .x, m = 1, n = n, n_days - 1,
                                                      cov_type = methods[i]))
    
    u_comb_mat[i,] = as.numeric(f_comb_agg - true_sep)
  }
  
  u_comb_mat
} 

aggForecast = function(f_vec, u_mat, out_idx, m, n, n_days, cov_type = "Sample"){
  u_in = t(u_mat[,-out_idx])
  mean_in = colMeans(u_in)
  
  if(cov_type == "EW"){
    return(mean(pull(f_vec)))
  }
  
  if(cov_type == "Sample"){
    cov_est = cov(u_in)
  }
  
  if(cov_type == "Linear"){
    cov_est = linshrink_cov(u_in)
  }
  
  if(cov_type == "Cor"){
    cov_est = getConstCorCov(u_in)
  }
  
  if(cov_type == "Var"){
    cov_est = diag(diag(cov(u_in)) + 1e-10) #some experts have 0 variance
  }
  
  if(cov_type == "S+EW"){
    w_SEW = SoptPlusEW(u_in)
    return(w_SEW %*% unlist(f_vec - mean_in))
  }
  
  if(cov_type == "Rob"){
    if(is_singular(cov(u_in))){ #to mitigate possible multicolinearity
      cov_est = diag(ncol(u_in))
    }else{
      cov_est = cov.rob(u_in,seed = 1:20240422)$cov + 1e-8 * diag(ncol(u_in)) #to prevent generate a singular covariance estimation
    }
  }
  
  if(m > 1){
    f_agg = t(getPoolWeights(solve(cov_est), m, n)) %*% unlist(f_vec - mean_in)
  }else{
    w = rowSums(solve(cov_est))
    w = w/sum(w)
    f_agg = w %*% unlist(f_vec - mean_in)
  }
  
  f_agg
}

getConstCorCov = function(u_in){
  n = nrow(u_in)
  k = ncol(u_in)
  
  sd_vec = sqrt(diag(cov(u_in)))
  cor_mat = cor(u_in)
  
  lower_tri_idx = which(lower.tri(cor_mat, diag = FALSE))
  
  rho_mean = mean(cor_mat[lower_tri_idx])
  
  rho_est = rho_mean/(1+2*k/(n-1) + 3/(n-1))
  
  const_cor_mat = diag(1-rho_est, k, k)  + rho_est
  
  cov_est = diag(sd_vec) %*% const_cor_mat %*% diag(sd_vec)
  
  return(cov_est)
}

SoptPlusEW <- function(train_data){
  X <- scale(as.matrix(train_data), center = TRUE, scale = FALSE)
  n = nrow(X) - 1
  k = ncol(X)
  ew <- rep(1/ncol(X), ncol(X)) %>% as.matrix()
  l <- matrix(rep(1, ncol(X)), ncol = 1)
  covar <- 1/n*t(X)%*%X + 1e-10 * diag(ncol(train_data))
  covar_inv <- solve(covar)
  ow <- as.numeric((covar_inv%*%l)/(t(l)%*%covar_inv%*%l %>% as.numeric()))
  V <- (n/(n-k+1))*((k-1)/(n-k))*(1/(t(l)%*%covar_inv%*%l))
  BB <- (t(l)%*%covar%*%l)/k^2-n/(n-k+1)*(1/(t(l)%*%covar_inv%*%l))
  B <- max(0, BB)
  lambda <- V/(B+V)
  weight <- c(lambda)*c(ew) + c(1-lambda)*ow
  return(weight)
}

getWRMSSE = function(u_mat, weight_vec = NULL){
  n_prods = nrow(u_mat)
  if(is.null(weight_vec)){
    weight_vec = rep(1/n_prods, n_prods)
  }
  
  sqrt(rowMeans(u_mat^2)) %*% weight_vec
}

extractNSub <- function(input_string) {
  # Use regular expression to extract the number after "Idx"
  match <- regmatches(input_string, regexpr("Idx(\\d+)", input_string))
  number <- as.numeric(sub("Idx", "", match))
  return(number)
}

