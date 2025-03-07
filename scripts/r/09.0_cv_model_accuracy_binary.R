################################################################################
################################################################################

# DESCRIPTION:

# This script summarizes classification performance across model types

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

# NOTES:

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpmisc)
library(viridis)
library(probably)
library(grid)
library(yardstick)

# 1.2 Parameters ---------------------------------------------------------------

ds_types = c('total', 'woody')
model_list = c('lgbm', 'rf', 'svmLinear', 'svmPoly', 'xgb')
metric_type = 'roc_auc' # Metric used to determine best hyperparameters
binary_threshold_metric = 'j_index' # Metric used to choose best threshold
threshold_deterimination_method = 'calculate' # Choose 'average' or 'calculate'

in_dir = paste0('output/08_cv_results/binary_loocv_10/')
out_dir = paste0('output/09_model_comparison/binary_loocv_10/accuracy_summary/')

# ==============================================================================
# 2. LOOP DATASETS AND MODELS AND CALCULATE OVERALL ACCURACY METRICS ===========
# ==============================================================================

mod_evaluation = data.frame()

for(ds_type in ds_types){ # START DATASET TYPE LOOP

  print(paste0('The dataset type is: ', ds_type))
  cat('\n')
  
  for(model_type in model_list){ # START MODEL TYPE LOOP
  
    print(paste0('The model type is: ', model_type))
    cat('\n')
    
    # 2.1 Gather and combine files ---------------------------------------------
    
    # Get lists of all files
    threshold_files = list.files(paste0(in_dir, '/thresholds/'), paste0('classification_thresholds_', ds_type, '_', model_type), full.names = TRUE)
    final_predictions_files = list.files(paste0(in_dir, '/predictions/'), paste0('final_predictions_', ds_type, '_', model_type), full.names = TRUE)

    # Read in as dataframes and combine
    thresholds = do.call(rbind, lapply(threshold_files, read.csv, header = TRUE))
    final_predictions = do.call(rbind, lapply(final_predictions_files, read.csv, header = TRUE))

    # 2.2 Choose presence threshold --------------------------------------------
    
    # Specify presence_absence response variable correctly
    final_predictions$presence_absence = factor(final_predictions$presence_absence, levels = c('presence', 'absence'))
    
    # Set up F measure metric with beta parameter
    f_meas_beta2 = metric_tweak("f_meas_beta2", f_meas, beta = 2)
    
    # Using predicted probabilities from each outer fold test set, calculate accuracy metrics anew for each potential presence threshold    
    thresholds = final_predictions %>%
      mutate(presence_absence = as.factor(presence_absence)) %>%
      threshold_perf(truth = presence_absence, 
                     estimate = prob_presence, 
                     thresholds = seq(0, 1, by = 0.05),
                     metrics = metric_set(yardstick::accuracy, yardstick::bal_accuracy, yardstick::kap, yardstick::f_meas, f_meas_beta2, yardstick::j_index, yardstick::mcc, yardstick::sens, yardstick::spec, yardstick::precision))
    names(thresholds) = gsub('[.]', '', names(thresholds))
    
    # Get best threshold from calculated results
    best_threshold = thresholds[thresholds$metric == binary_threshold_metric,] %>%
      group_by(threshold) %>%
      summarise(estimate_avg = mean(estimate)) %>%
      slice_max(estimate_avg, n = 1) %>% # Choose threshold with maximum metric (j_index, sensitivity or specificity)
      slice_min(threshold, n = 1) # If there is a tie, choose the lower threshold
    
    # 2.3 Calculate accuracy metrics -------------------------------------------
    
    # Assign factor
    final_predictions$presence_absence = as.factor(final_predictions$presence_absence)
    final_predictions$presence_absence = factor(final_predictions$presence_absence, levels = c('presence', 'absence'))
    
    # Assign predictions
    final_predictions = final_predictions %>%
      mutate(
        predicted = make_two_class_pred(
          levels = levels(presence_absence), # A character vector of class levels
          estimate = prob_presence, # A single numeric vector corresponding to the class probabilities of the first level in levels
          threshold = best_threshold$threshold
        )
      )
    
    # Assign factor
    final_predictions$predicted = as.factor(final_predictions$predicted)
    
    # Calculate accuracy metrics
    fit_accuracy = round(yardstick::accuracy_vec(final_predictions$presence_absence, final_predictions$predicted), 2)
    fit_j_index = round(yardstick::j_index_vec(final_predictions$presence_absence, final_predictions$predicted), 2)
    fit_f_meas = round(yardstick::f_meas_vec(final_predictions$presence_absence, final_predictions$predicted), 2)
    fit_kap = round(yardstick::kap_vec(final_predictions$presence_absence, final_predictions$predicted), 2)
    
    # 2.4 Tidy and combine -----------------------------------------------------
  
    mod_df = data.frame(ds_type = ds_type, 
                        mod = model_type, 
                        threshold = best_threshold$threshold,
                        accuracy = fit_accuracy,
                        j_index = fit_j_index,
                        f_meas = fit_f_meas,
                        kappa = fit_kap)
    
    mod_evaluation = rbind(mod_evaluation, mod_df)
    
  } # END MODEL TYPE LOOP
  
} # END DATASET TYPE LOOP

# ==============================================================================
# 3. SAVE ======================================================================
# ==============================================================================

write.csv(mod_evaluation, paste0(out_dir, 'binary_model_accuracy_comparison.csv'), row.names = FALSE)
