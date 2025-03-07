################################################################################
################################################################################

# DESCRIPTION:

# This script tests different thresholds for assigning biomass 'presence'
# based on presence proability results and chooses best threshold

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(probably)

source('scripts/r/000.0_yardstick_bias.R')
source('scripts/r/000.1_yardstick_rmspe.R')

# 1.2 Parameters ---------------------------------------------------------------

# R
in_version = 'v20240514'
out_version = 'v20240514'
ds_types = c('total', 'woody')
ds_type = 'woody'
mod_name = 'rf'
binary_threshold_metric = 'j_index'

# Directories
in_dir = 'output/11_mc_final_models/binary/'
out_dir = 'output/11_mc_final_models/binary/thresholds/'

# 1.3 Report -------------------------------------------------------------------

cat('\n')
print(paste0('In version: ', in_dir))
print(paste0('Out version: ', out_dir))
print(paste0('Data set type: ', ds_type))
print(paste0('Model: ', mod_name))
print(paste0('The presence probability threshold metric is: ', binary_threshold_metric))
cat('\n')

################################################################################
################################################################################

# ==============================================================================
# 2. TEST THRESHOLDS ===========================================================
# ==============================================================================

for(ds_type in ds_types){ # START DATASET TYPE LOOP
  
  print(paste0('The dataset type is: ', ds_type))
  cat('\n')
  
  # Set up data frame
  thresholds_df = data.frame()
  best_thresholds_df = data.frame()
  
  for(i in 1:100){ # START MONTE CARLO ITERATION LOOP
  
    print(paste0('The Monte Carlo iteration is: ', i))
    cat('\n')
    
    # 2.1 Read in model predictions ----------------------------------------------
    
    # Read in model predictions
    predictions_df = read.csv(paste0(in_dir, 'final_predictions_', ds_type, '_binary_', out_version, '_', i, '.csv'))
    
    # Tidy
    predictions_df$presence_absence = factor(predictions_df$presence_absence, levels = c('presence', 'absence'))
    
    # 2.2 Test different thresholds ----------------------------------------------
  
    # Set up F measure metric with beta parameter
    f_meas_beta2 = metric_tweak("f_meas_beta2", f_meas, beta = 2)
    
    # Calculate threshold performances
    thresholds = predictions_df %>%
      threshold_perf(truth = presence_absence, 
                     estimate = prob_presence, 
                     thresholds = seq(0, 1, by = 0.05),
                     metrics = metric_set(yardstick::accuracy, yardstick::bal_accuracy, yardstick::kap, yardstick::f_meas, f_meas_beta2, yardstick::j_index, yardstick::mcc, yardstick::sens, yardstick::spec, yardstick::precision))
    
    # Tidy threshold data
    thresholds = thresholds %>% rename(threshold = .threshold, metric = .metric, estimator = .estimator, estimate = .estimate)
  
    # Append
    thresholds$mod = mod_name
    thresholds$ds_type = ds_type
    thresholds$mc = i
    thresholds_df = bind_rows(thresholds_df, thresholds)
    
    # 2.3 Choose best threshold --------------------------------------------------
    
    # Find best threshold
    best_threshold = thresholds[thresholds$metric == binary_threshold_metric,] %>%
      group_by(threshold) %>%
      summarise(estimate_avg = mean(estimate)) %>%
      slice_max(estimate_avg, n = 1) %>% # Choose threshold with maximum metric (j_index, sensitivity or specificity)
      slice_min(threshold, n = 1) # If there is a tie, choose the lower threshold
    
    # Append
    best_thresholds_df = bind_rows(best_thresholds_df, data.frame(mod = mod_name, ds_type = ds_type, mc = i, threshold = best_threshold$threshold))
  
  } # END MONTE CARLO ITERATION LOOP
  
  write.csv(thresholds_df, paste0(out_dir, 'presence_thresholds_', ds_type, '_binary_', out_version, '.csv'), row.names = FALSE)
  write.csv(best_thresholds_df, paste0(out_dir, 'presence_thresholds_',  binary_threshold_metric, '_', ds_type, '_binary_', out_version, '.csv'), row.names = FALSE)

} # END DATASET TYPE LOOP