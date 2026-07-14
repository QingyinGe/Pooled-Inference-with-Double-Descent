## ============================================================================
## GoodnessOfFit.R
##
## Tests whether the log-variance of pairwise expert forecast-error contrasts
## decomposes additively into a "product" effect and an "expert-contrast"
## effect, via a two-way-fixed-effects (TWFE) regression:
##
##     log( var(x_{n,m}) ) = grand_mean + row_effect[n] + col_effect[m] + resid
##
## where m indexes "products" (an M4 time series, or an HHS region for flu)
## and n indexes one of the (n_experts - 1) orthogonal contrasts among the
## experts forecasting that product (built from Gamma0, see getGamma0() in
## func.R). The R^2 / adjusted R^2 of that regression is the "goodness of
## fit" statistic reported below.
##
## This single file supports two datasets:
##   - "flu": CDC flu-forecast data, one adjusted R^2 PER SEASON (8 seasons)
##            -> bar chart
##   - "m4" : M4 competition data, one adjusted R^2 PER SCENARIO (68 monthly
##            + 2 daily = 70 scenarios) -> boxplot on a log10 y-axis
##
## Set `dataset_type` below and run the whole file; the matching pipeline
## loads its data, computes the R^2 values, prints/saves a table, and draws
## + saves its plot.
## ============================================================================

library(tidyverse)
library(Matrix)
library(data.table)

# dataset_type controls which pipeline runs at the bottom of this file.
dataset_type = "flu"   # "flu" or "m4"

# getGamma0() lives in func.R (assumes this script is run with the working
# directory set to Flu-forecast/, where a copy of func.R lives alongside it).
source("func.R")

theme_slides = theme(text = element_text(size = 15), legend.position = "top")


## ----------------------------------------------------------------------------
## SHARED CORE (identical math for both datasets)
## ----------------------------------------------------------------------------

## Closed-form TWFE R^2 (and adjusted R^2) for a *complete* balanced two-way
## grid of variances, var_mat[n, m] (n = 1..n_contrasts contrasts, m = 1..
## n_prods products). Because the grid is a complete rectangle (every n
## crossed with every m, no missing cells), this closed form is algebraically
## identical to lm(log(var_mat) ~ factor(m) + factor(n)).
getTWFE_R2 = function(var_mat){
  log_mat = log(var_mat)
  n_contrasts = nrow(var_mat)
  n_prods = ncol(var_mat)

  row_means = rowMeans(log_mat)    # per contrast n
  col_means = colMeans(log_mat)    # per product m
  grand_mean = mean(log_mat)

  fitted = outer(row_means, col_means, function(r, c) r + c - grand_mean)
  resid = log_mat - fitted

  ss_tot = sum((log_mat - grand_mean)^2)
  ss_res = sum(resid^2)
  r2 = 1 - ss_res / ss_tot

  # n_obs = number of (contrast, product) cells; params = dummy variables for
  # each factor with one reference level dropped per factor (matches lm()'s
  # default dummy-coding), plus 1 for the intercept.
  n_obs = n_prods * n_contrasts
  n_params = (n_prods - 1) + (n_contrasts - 1)
  df_resid = n_obs - n_params - 1
  adj_r2 = 1 - (1 - r2) * (n_obs - 1) / df_resid

  c(R2 = r2, adjR2 = adj_r2)
}

## Turn a raw forecast-error matrix into the n_contrasts x n_prods variance
## matrix that getTWFE_R2() needs.
##
## err_mat must be laid out as n_prods stacked blocks of n_experts rows each
## (block m = the n_experts forecast-error time series for product m, in a
## fixed expert order), with n_periods columns (time). Row order within each
## block must be consistent across products since Gamma0 is applied per block.
computeVarMat = function(err_mat, n_experts, n_prods){
  n_periods = ncol(err_mat)

  Gamma0 = getGamma0(n_experts)                             # n_experts x (n_experts-1) contrast basis
  Gamma = Matrix(diag(n_prods), sparse = TRUE) %x% Gamma0   # block-diagonal, one Gamma0 per product

  # Project each product's n_experts-dim error vector onto the
  # (n_experts-1)-dim contrast space, for every time period at once.
  X_mat = as.matrix(t(err_mat) %*% Gamma)              # n_periods x (n_prods * n_contrasts)
  X_mat = scale(X_mat, center = TRUE, scale = FALSE)   # demean each contrast-product column over time

  var_vec = colSums(X_mat^2) / n_periods               # var of each contrast-product column
  matrix(var_vec, nrow = n_experts - 1, ncol = n_prods) # reshape: rows = contrasts, cols = products
}


## ==============================================================================
## FLU PIPELINE
## One adjusted R^2 per flu season (8 seasons) -> bar chart.
## ==============================================================================
runFluGoodnessOfFit = function(csv_path = "point_ests_adj-w20172018.csv"){

  # Keep 1-week-ahead forecasts and the 10 HHS regions (drop the "US
  # National" aggregate). Exclude ReichLab_kde, UTAustin_edm (has a 13-week
  # submission gap in the 2017/2018 season), and the six ensemble/weighting
  # pseudo-models, so every remaining model has a complete panel every season.
  df = fread(csv_path) %>%
    filter(target == "1 wk ahead", location != "US National",
           !model_name %in% c("ReichLab_kde", "UTAustin_edm",
                               "constant-weights", "equal-weights",
                               "target-and-region-based-weights",
                               "target-based-weights",
                               "target-type-based-weights"))

  locations = sort(unique(df$location))
  n_prods = length(locations)          # 10 HHS regions ("products")
  seasons = sort(unique(df$Season))    # 8 flu seasons

  # Pivot one season's forecast errors ("err") to a wide (location, model)
  # x week matrix, rows ordered location-major / model_name-minor so each
  # location's block of n_experts rows is contiguous, as computeVarMat() needs.
  prepareFluSeason = function(season){
    df %>%
      filter(Season == season) %>%
      dplyr::select(location, model_name, Model.Week, err) %>%
      arrange(Model.Week) %>%
      pivot_wider(names_from = Model.Week, values_from = err) %>%
      arrange(location, model_name)
  }

  results = list()
  for(season in seasons){
    wide = prepareFluSeason(season)
    n_experts = length(unique(wide$model_name))

    err_mat = wide %>% dplyr::select(-location, -model_name) %>% as.matrix()
    err_mat = err_mat - rowMeans(err_mat)   # demean each (location, model) series over time first

    var_mat = computeVarMat(err_mat, n_experts, n_prods)
    r2_vec = getTWFE_R2(var_mat)

    results[[season]] = data.frame(
      Season = season, n_experts = n_experts, n_prods = n_prods,
      n_periods = ncol(err_mat), R2 = r2_vec["R2"], adjR2 = r2_vec["adjR2"]
    )

    cat(sprintf("Flu season %s: n_experts=%d n_prods=%d weeks=%d R2=%.4f adjR2=%.4f\n",
                season, n_experts, n_prods, ncol(err_mat), r2_vec["R2"], r2_vec["adjR2"]))
  }

  r2_df = do.call(rbind, results)
  rownames(r2_df) = NULL
  r2_df$Season = factor(r2_df$Season, levels = seasons)

  save(r2_df, file = "FluGoodnessOfFit_bySeason.RData")
  write.csv(r2_df, "FluGoodnessOfFit_bySeason.csv", row.names = FALSE)

  # Bar chart: one adjusted R^2 per season (a boxplot doesn't apply here
  # since there's only a single number per season, not a distribution).
  p_bar = ggplot(r2_df, aes(x = Season, y = adjR2)) +
    geom_col(fill = "grey70", color = "black") +
    geom_text(aes(label = round(adjR2, 3)), vjust = -0.5, size = 4) +
    theme_bw() + theme_slides +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1)) +
    ylab("Adjusted R2") + xlab("Season")

  print(p_bar)
  ggsave("FluGoodnessOfFit_bySeason_bar.pdf", plot = p_bar, width = 9, height = 6)

  r2_df
}


## ==============================================================================
## M4 PIPELINE
## One adjusted R^2 per scenario (68 monthly starting-date clusters + 2 daily
## clusters = 70 scenarios) -> boxplot on a log10 y-axis.
## ==============================================================================
runM4GoodnessOfFit = function(n_experts = 17, scaled_data_dir = "../Pooled-Inference-with-Double-Descent-DC18/ScaledData"){

  # A "scenario" = the group of M4 series sharing the rank_idx-th most common
  # starting date within a given frequency (Monthly/Daily); each such series
  # is one "product". This mirrors the original Code/GoodnessOfFit.R logic.
  getTimeSeriesNames = function(data_freq, rank_idx){
    info = read.csv(file.path(scaled_data_dir, "M4-info.csv"))
    info_sub = info %>% filter(SP == data_freq)
    top_dates = names(sort(table(info_sub$StartingDate), decreasing = TRUE))
    chosen = info_sub %>% filter(StartingDate == top_dates[rank_idx])
    chosen$M4id
  }

  # Build the block-stacked error matrix for one scenario: n_prods products
  # (M4 series), n_experts consecutive rows per product, n_periods
  # forecast-horizon columns (F2, F3, ... in the err file).
  prepareM4Scenario = function(data_freq, rank_idx, u_data){
    ids = getTimeSeriesNames(data_freq, rank_idx)
    n_prods = length(ids)
    err_mat = as.matrix(
      u_data %>% filter(id %in% ids) %>%
        arrange(as.numeric(gsub("\\D", "", id))) %>%
        dplyr::select(-id, -group_id)
    )
    list(err_mat = err_mat, n_prods = n_prods)
  }

  results = list()

  for(data_freq in c("Monthly", "Daily")){
    u_data = fread(file.path(scaled_data_dir, paste0(data_freq, "_err.csv"))) %>% dplyr::select(-V1)
    n_scenarios = if(data_freq == "Monthly") 68 else 2

    for(rank_idx in 1:n_scenarios){
      prep = prepareM4Scenario(data_freq, rank_idx, u_data)
      var_mat = computeVarMat(prep$err_mat, n_experts, prep$n_prods)
      r2_vec = getTWFE_R2(var_mat)

      scenario = paste0(data_freq, rank_idx)
      results[[scenario]] = data.frame(
        scenario = scenario, data_freq = data_freq, rank_idx = rank_idx,
        n_prods = prep$n_prods, n_periods = ncol(prep$err_mat),
        R2 = r2_vec["R2"], adjR2 = r2_vec["adjR2"]
      )

      cat(sprintf("%s: n_prods=%d R2=%.4f adjR2=%.4f\n",
                  scenario, prep$n_prods, r2_vec["R2"], r2_vec["adjR2"]))
    }
  }

  r2_df = do.call(rbind, results)
  rownames(r2_df) = NULL

  save(r2_df, file = "M4GoodnessOfFit_70scenarios.RData")
  write.csv(r2_df, "M4GoodnessOfFit_70scenarios.csv", row.names = FALSE)

  # Boxplot of the ~70 adjusted R^2 values (linear y-axis). Saved as PDF.
  p_box = ggplot(r2_df, aes(x = "", y = adjR2)) +
    geom_boxplot(outlier.alpha = 0.4, width = 0.3) +
    theme_bw() + theme_slides +
    xlab(NULL) + ylab("Adjusted R2")

  print(p_box)
  ggsave("M4GoodnessOfFit_boxplot.pdf", plot = p_box, width = 8, height = 6)

  r2_df
}


## ----------------------------------------------------------------------------
## RUN
## ----------------------------------------------------------------------------
if(dataset_type == "flu"){
  r2_df = runFluGoodnessOfFit()
} else if(dataset_type == "m4"){
  r2_df = runM4GoodnessOfFit()
} else {
  stop('dataset_type must be "flu" or "m4"')
}
