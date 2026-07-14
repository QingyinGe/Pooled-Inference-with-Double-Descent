library(tidyverse)
library(nlshrink)
library(MASS)
library(Matrix)
library(data.table)
# library(pls)
source("Code/func.R")

# Delete the ones with missing data

target_week = "1 wk ahead"; PCA_threshold = 0.8 # Figure ?
# target_week = "1 wk ahead"; PCA_threshold = 0.7 # Figure ?
# target_week = "2 wk ahead"; PCA_threshold = 0.8 # Figure ?

df = fread("CleanedData/FluForecasting/point_ests_adj-w20172018.csv") %>% 
  filter(target == target_week & location != "US National" & 
           !model_name %in% c("ReichLab_kde", "UTAustin_edm", "constant-weights",
                              "equal-weights", "target-and-region-based-weights", 
                              "target-based-weights", "target-type-based-weights"))

# sum(is.na(df))

n_experts = length(unique(df$model_name))
n_variables = length(unique(df$location))
n_year = 7
methods = c("S:Linear", "S:Cor", "S:Var")
n_methods = length(methods)

ratio_mat = NULL

for(i in 1:(n_year-3)){
  df_train = cbind(prepareFluData(df, paste0(2009 + i, '/', 2010 + i)),
                   prepareFluData(df, paste0(2010 + i, '/', 2011 + i)),
                   prepareFluData(df, paste0(2011 + i, '/', 2012 + i)),
                   prepareFluData(df, paste0(2012 + i, '/', 2013 + i)))
  
  df_test = prepareFluData(df, paste0(2013 + i, '/', 2014 + i))
  
  # stopifnot(df_train$location == df_test$location, 
  #           df_train$model_name == df_test$model_name)
  
  df_train = df_train %>% dplyr::select(-location, -model_name)
  df_test = df_test %>% dplyr::select(-location, -model_name)
  
  Y_train = t(apply(df_train, 2, function(x) sumEachProd(x, n_variables, n_experts)/n_experts))
  Y_test = t(apply(df_test, 2, function(x) sumEachProd(x, n_variables, n_experts)/n_experts))
  
  Gamma0 = getGamma0(n_experts)
  Gamma = Matrix(diag(n_variables), sparse = TRUE) %x% Gamma0
  
  X_train = as.matrix(t(df_train) %*% Gamma)
  X_test = as.matrix(t(df_test) %*% Gamma)
  
  PCR_test_err = getPCRTestErr(X_train, Y_train, X_test, Y_test, threshold = PCA_threshold)
  ridge0_test_err_std = getRidge0TestErr(X_train, Y_train, X_test, Y_test, stdize = TRUE)
  
  pool_idx = rep(1:n_experts, n_variables) + 
    rep(0:(n_variables - 1) * n_experts, each = n_experts)
  
  df_test_debiased = df_test - rowMeans(df_train)
  
  u_pool_Linear = getFastPoolLinearTestU(t(df_train), t(df_test_debiased), n_variables,
                                         n_experts, n_obs = ncol(df_train))
  
  L2_pool = apply(u_pool_Linear, 1, function(x) mean(x^2))
  # L2_pool = NA
  L2_PCR = apply(PCR_test_err, 1, function(x) mean(x^2))
  L2_ridge0_std = apply(ridge0_test_err_std, 1, function(x) mean(x^2))
  L2_EW = apply(Y_test, 2, function(x) mean(x^2))
  
  
  # Separate methods --------------------------------------------------------
  
  L2_sep_methods = matrix(NA, n_methods, n_variables, dimnames = list(methods))
  
  for(j in 1:n_variables){
    sep_idx = 1:n_experts + (j-1) * n_experts
    
    df_train_sep = t(df_train[sep_idx,])
    df_test_sep = t(df_test[sep_idx,])
    
    L2_sep_methods[,j] = map_dbl(1:n_methods, ~ getSepTestL2(df_train_sep, df_test_sep, methods[.x]))
  }
  
  
  # print(round(cbind(L2_ridge0_std, L2_pool, L2_EW, t(L2_sep_methods),L2_PCR),2))
  
  ratio_mat = rbind(ratio_mat,
    cbind(L2_ridge0_std, L2_pool, L2_EW, t(L2_sep_methods))/L2_PCR)
}


colnames(ratio_mat) = c("P:Ridgeless", "P:Linear", "S:EW", methods)

summary(ratio_mat)

theme_slides = theme(text=element_text(size=15),
                     legend.position = "top")

df4plot_long = pivot_longer(as.data.frame(ratio_mat), cols = `P:Ridgeless`:`S:Var`,
                            values_to = "Ratio", names_to = "Method")
df4plot_long$Method = factor(df4plot_long$Method, 
                             levels = c("P:Ridgeless", "P:Linear", "S:EW", "S:Var", "S:Cor", "S:Linear"))

plot1 = ggplot(df4plot_long %>% filter(!(Method %in% c("P:Ridgeless", "P:Linear"))), aes(x = Method, y = Ratio)) +
  geom_boxplot(outlier.alpha = 0.05) + coord_cartesian(ylim = c(0.1,10)) + 
  geom_hline(yintercept = 1, col = 'red') + 
  theme_bw() + theme_slides + ylab("Loss Ratio (Benchmarks/P:PCR)") + xlab("Benchmarks") + 
  scale_y_log10() + 
  annotation_logticks(sides = "l")

print(plot1)

if(target_week == "1 wk ahead" & PCA_threshold == 0.8){
  pdf("Figures/flu1wk80.pdf")
}else if(target_week == "1 wk ahead" & PCA_threshold == 0.7){
  pdf("Figures/flu1wk70.pdf")
}else if (target_week == "2 wk ahead" & PCA_threshold == 0.8){
  pdf("Figures/flu2wk80.pdf")
}
print(plot1)
dev.off()

