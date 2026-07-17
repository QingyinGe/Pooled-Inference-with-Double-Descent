## ============================================================================
## Tests whether the log-variance of pairwise expert forecast-error contrasts
## decomposes additively into a "product" effect and an "expert-contrast"
## effect, via a two-way-fixed-effects (TWFE) regression:
##
##     log( var(x_{n,m}) ) = grand_mean + row_effect[n] + col_effect[m] + resid
##
## Runs both datasets in turn:
##   - flu: CDC flu-forecast data, one adjusted R^2 PER SEASON (8 seasons),
##          computed separately for the 1 wk ahead and 2 wk ahead targets
##          -> bar chart (1 wk ahead) + printed 1wk-vs-2wk table
##   - M4 : M4 competition data, one adjusted R^2 PER SCENARIO (68 monthly
##          + 2 daily = 70 scenarios) -> histogram
## ============================================================================

library(tidyverse)
library(Matrix)
library(data.table)

source("Code/func.R")
theme_slides = theme(text = element_text(size = 15), legend.position = "top")

## ----------------------------------------------------------------------------
## Flu dataset
## ----------------------------------------------------------------------------

csv_path = "CleanedData/FluForecasting/point_ests_adj-w20172018.csv"

# target_week = "1 wk ahead"
target_week = "2 wk ahead" 

df = fread(csv_path) %>%
  filter(target == target_week, location != "US National",
         !model_name %in% c("ReichLab_kde", "UTAustin_edm",
                             "constant-weights", "equal-weights",
                             "target-and-region-based-weights",
                             "target-based-weights",
                             "target-type-based-weights"))

locations = sort(unique(df$location))
n_prods = length(locations)          # 10 HHS regions ("products")
seasons = sort(unique(df$Season))    # 8 flu seasons

results = list()
for(season in seasons){
  wide = prepareFluData(df, season)
  n_experts = length(unique(wide$model_name))

  err_mat = wide %>% dplyr::select(-location, -model_name) %>% as.matrix()

  var_mat = computeVarMat(err_mat, n_experts, n_prods)
  # getTWFE_R2() fits log(var_mat) ~ region + contrast and returns that
  # fit's R^2 and adjusted R^2.
  r2_vec = getTWFE_R2(var_mat)

  results[[season]] = data.frame(
    Season = season, n_experts = n_experts, n_prods = n_prods,
    n_periods = ncol(err_mat), R2 = r2_vec["R2"], adjR2 = r2_vec["adjR2"]
  )

  cat(sprintf("Flu %s season %s: n_experts=%d n_prods=%d weeks=%d R2=%.4f adjR2=%.4f\n",
              target_week, season, n_experts, n_prods, ncol(err_mat), r2_vec["R2"], r2_vec["adjR2"]))
}

# r2_df_flu_1wk = do.call(rbind, results)  ## Table EC.2 second column
# rownames(r2_df_flu_1wk) = NULL
# r2_df_flu_1wk$Season = factor(r2_df_flu_1wk$Season, levels = seasons)

r2_df_flu_2wk = do.call(rbind, results)  ## Table EC.2 third column
rownames(r2_df_flu_2wk) = NULL
r2_df_flu_2wk$Season = factor(r2_df_flu_2wk$Season, levels = seasons)


## ----------------------------------------------------------------------------
## M4 dataset
## ----------------------------------------------------------------------------

n_experts = 17
scaled_data_dir = "CleanedData/M4ScaledData"

results = list()

for(data_freq in c("Monthly", "Daily")){
  u_data = fread(file.path(scaled_data_dir, paste0(data_freq, "_err.csv"))) %>% dplyr::select(-V1)
  n_scenarios = if(data_freq == "Monthly") 68 else 2

  for(rank_idx in 1:n_scenarios){
    # prepareM4Scenario() builds the block-stacked error matrix for one
    # scenario: n_prods M4 series, n_experts consecutive rows per series.
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

r2_df_m4 = do.call(rbind, results)
rownames(r2_df_m4) = NULL

save(r2_df_m4, file = "Result/M4GoodnessOfFit_70scenarios.RData")
write.csv(r2_df_m4, "Result/M4GoodnessOfFit_70scenarios.csv", row.names = FALSE)


# Histogram of the ~70 adjusted R^2 values (linear x-axis). Saved as PDF.
p_hist = ggplot(r2_df_m4, aes(x = adjR2)) +
  geom_histogram(bins = 15, color = "black", fill = "white") +
  theme_bw() + theme_slides +
  xlab("Adjusted R2") + ylab("Count")

print(p_hist)
ggsave("Result/M4GoodnessOfFit_histogram.pdf", plot = p_hist, width = 8, height = 6)
