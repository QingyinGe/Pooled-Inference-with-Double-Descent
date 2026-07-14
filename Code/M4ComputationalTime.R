# This code generates Table 3. --------------------------------------------

library(tidyverse)
library(Matrix)
library(MCMCpack)
library(data.table)
library(doParallel)
library(doRNG)
library(xtable)
pkg_names = c("tidyverse", "nlshrink", "MASS")

source("Code/func.R")

time_mat = matrix(NA, nrow = 8, ncol = 2, 
                  dimnames = list(c(paste0("Monthly", 1:7), "Daily1"), c("P:Linear", "P:Ridgeless")))

# Monthly Computational Time ----------------------------------------------
data_freq = "Monthly"

f_data = fread(paste0("CleanedData/M4ScaledData/", data_freq, "_forecast.csv")) %>% dplyr::select(-V1)
true_data = fread(paste0("CleanedData/M4ScaledData/", data_freq, "_truth.csv")) %>% dplyr::select(-V1)
u_data = fread(paste0("CleanedData/M4ScaledData/", data_freq, "_err.csv")) %>% dplyr::select(-V1)

rank_vec = 1:7

for(i in rank_vec){
  time_series_names = getTimeSeriesNames(data_freq, i)
  
  n_prods = length(time_series_names)
  # n_prods = 100
  n_experts = 17
  # n_periods = 18
  n_periods = ncol(true_data) - 1
  
  true_comb = as.data.frame(true_data %>% filter(id %in% time_series_names[1:n_prods]) %>% 
                              arrange(as.numeric(gsub("\\D", "", id))) %>% dplyr::select(-id))
  f_comb = as.data.frame(f_data %>% filter(id %in% time_series_names[1:n_prods]) %>% 
                           arrange(as.numeric(gsub("\\D", "", id))) %>% dplyr::select(-id, -group_id))
  u_comb = as.data.frame(u_data %>% filter(id %in% time_series_names[1:n_prods]) %>% 
                           arrange(as.numeric(gsub("\\D", "", id))) %>% dplyr::select(-id, -group_id))
  
  # P:Linear ----------------------------------------------------------------
  pool_idx = rep(1:n_experts, n_prods) + rep(0:(n_prods - 1) * n_experts, each = n_experts)
  
  # The following is too slow when n_prods is large
  # u_pool_Linear1 = getPooledU(pool_idx, "Linear", n_prods, n_experts, f_comb, u_comb, true_data)
  
  start_time = Sys.time()
  u_pool_Linear = getPooledUFast(pool_idx, n_prods, n_experts, u_comb, true_comb)
  end_time = Sys.time()
  
  time_mat[i,1] = as.numeric(difftime(end_time, start_time, units = "secs"))/ncol(true_comb)
  
  L2_pool = apply(u_pool_Linear, 1, function(x) mean(x^2))
  
  
  # P:Ridge0+ ---------------------------------------------------------------
  start_time = Sys.time()
  Y_mat = t(apply(u_comb, 2, function(x) sumEachProd(x, n_prods, n_experts)/n_experts))
  
  Gamma0 = getGamma0(n_experts)
  Gamma = Matrix(diag(n_prods), sparse = TRUE) %x% Gamma0
  
  X_mat = as.matrix(t(u_comb) %*% Gamma)
  
  ridge0_test_err_std = map_dfc(1:nrow(X_mat), ~ getRidge0TestErr(X_mat[-.x,], Y_mat[-.x,], X_mat[.x,,drop = FALSE],
                                                                  Y_mat[.x,,drop = FALSE], stdize = TRUE))
  end_time = Sys.time()
  
  time_mat[i,2] = as.numeric(difftime(end_time, start_time, units = "secs"))/ncol(true_comb)
  
  L2_ridge0_std = apply(ridge0_test_err_std, 1, function(x) mean(x^2))
}



# Daily Computational Time ------------------------------------------------
data_freq = "Daily"

f_data = fread(paste0("CleanedData/M4ScaledData/", data_freq, "_forecast.csv")) %>% dplyr::select(-V1)
true_data = fread(paste0("CleanedData/M4ScaledData/", data_freq, "_truth.csv")) %>% dplyr::select(-V1)
u_data = fread(paste0("CleanedData/M4ScaledData/", data_freq, "_err.csv")) %>% dplyr::select(-V1)


time_series_names = getTimeSeriesNames(data_freq, 1)

n_prods = length(time_series_names)
# n_prods = 100
n_experts = 17
# n_periods = 18
n_periods = ncol(true_data) - 1

true_comb = as.data.frame(true_data %>% filter(id %in% time_series_names[1:n_prods]) %>% 
                            arrange(as.numeric(gsub("\\D", "", id))) %>% dplyr::select(-id))
f_comb = as.data.frame(f_data %>% filter(id %in% time_series_names[1:n_prods]) %>% 
                         arrange(as.numeric(gsub("\\D", "", id))) %>% dplyr::select(-id, -group_id))
u_comb = as.data.frame(u_data %>% filter(id %in% time_series_names[1:n_prods]) %>% 
                         arrange(as.numeric(gsub("\\D", "", id))) %>% dplyr::select(-id, -group_id))

# P:Linear ----------------------------------------------------------------
pool_idx = rep(1:n_experts, n_prods) + rep(0:(n_prods - 1) * n_experts, each = n_experts)

# The following is too slow when n_prods is large
# u_pool_Linear1 = getPooledU(pool_idx, "Linear", n_prods, n_experts, f_comb, u_comb, true_data)

start_time = Sys.time()
u_pool_Linear = getPooledUFast(pool_idx, n_prods, n_experts, u_comb, true_comb)
end_time = Sys.time()

time_mat[8,1] = as.numeric(difftime(end_time, start_time, units = "secs"))/ncol(true_comb)

L2_pool = apply(u_pool_Linear, 1, function(x) mean(x^2))


# P:Ridge0+ ---------------------------------------------------------------
start_time = Sys.time()
Y_mat = t(apply(u_comb, 2, function(x) sumEachProd(x, n_prods, n_experts)/n_experts))

Gamma0 = getGamma0(n_experts)
Gamma = Matrix(diag(n_prods), sparse = TRUE) %x% Gamma0

X_mat = as.matrix(t(u_comb) %*% Gamma)

ridge0_test_err_std = map_dfc(1:nrow(X_mat), ~ getRidge0TestErr(X_mat[-.x,], Y_mat[-.x,], X_mat[.x,,drop = FALSE],
                                                                Y_mat[.x,,drop = FALSE], stdize = TRUE))
end_time = Sys.time()

time_mat[8,2] = as.numeric(difftime(end_time, start_time, units = "secs"))/ncol(true_comb)

L2_ridge0_std = apply(ridge0_test_err_std, 1, function(x) mean(x^2))


# Print out results -------------------------------------------------------
time_mat = time_mat[c(paste0("Monthly", 7:5), "Daily1", paste0("Monthly", 4:1)),]

time_mat = cbind(NPooled = c(1048, 1225, 1285, 1506, 1517, 2643, 6258, 6782),
                 time_mat)

xtable(t(time_mat), digit = 2)
