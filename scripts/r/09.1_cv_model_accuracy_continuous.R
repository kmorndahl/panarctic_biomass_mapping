################################################################################
################################################################################

# DESCRIPTION:

# This script summarizes regression performance across model types

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
library(yardstick)
library(classInt)

source('scripts/000.0_yardstick_bias.R')
source('scripts/000.1_yardstick_rmspe.R')

# 1.2 Parameters ---------------------------------------------------------------

ds_types = c('total', 'woody')
model_list = c('linearBaseline_log', 'lgbm', 'lgbm_log', 'rf', 'rf_log', 'svmLinear', 'svmLinear_log', 'svmPoly', 'svmPoly_log', 'xgb', 'xgb_log')
ds_subsets = c('all', 'presence')
metric_type = 'rmse' # Metric used to determine best hyperparameters

in_dir = paste0('output/08_cv_results/continuous_loocv_10/')
out_dir = paste0('output/09_model_comparison/continuous_loocv_10/accuracy_summary/')

################################################################################
################################################################################

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
    
    train_predictions_files = list.files(paste0(in_dir, '/predictions/'), paste0('final_predictions_train_', ds_type, '_', model_type), full.names = TRUE)
    final_predictions_files = list.files(paste0(in_dir, '/predictions_train/'), paste0('final_predictions_', ds_type, '_', model_type), full.names = TRUE)

    # Remove log transformed models if not specified
    if(!grepl('log', model_type)){
      train_predictions_files = train_predictions_files[!grepl('_log', train_predictions_files)]
      final_predictions_files = final_predictions_files[!grepl('_log', final_predictions_files)]
    }
    
    train_predictions = do.call(rbind, lapply(train_predictions_files, read.csv, header = TRUE))
    final_predictions = do.call(rbind, lapply(final_predictions_files, read.csv, header = TRUE))

    print('Data read in')
    
    for(ds_subset in ds_subsets){ # START ZERO HANDLING LOOP
      
      print(paste0('The data subset type is: ', ds_subset))
      cat('\n')
      
      # 2.2 Subset data if specified -------------------------------------------
      
      # Remove user identified non-vegetated observations if presence only subset specified
      if(ds_subset != 'all'){
        final_predictions_subset = final_predictions[final_predictions$dataset_id != 'nonveg',]
      }else{
        final_predictions_subset = final_predictions
      }
      
      print('Data subset (if specified)')
      
      # 2.3 Apply bias correction to log transformed models  -------------------
      
      if(grepl('log', model_type)){
          
          # Calculate correction factor - binned smear
          train_predictions = train_predictions %>%
            group_by(outer_fold) %>%
            mutate(response_bin_breaks = list(c(head(log(classIntervals(exp(predicted), n=5, style="fisher", warnLargeN = FALSE)$brks), -1), Inf))) %>% # Calculates breaks, replaces last value with Inf to accommodate values in test predictions that exceed those in training predictions
            mutate(response_bin = cut(predicted, breaks = response_bin_breaks[[1]], include.lowest = TRUE)) # Assigns bins based on breaks
          train_predictions = train_predictions %>%
            group_by(outer_fold, response_bin) %>% # Group by outer fold AND response bin
            mutate(correction_factor = sum(exp(biomass_density_gm2  - predicted))/length(biomass_density_gm2)) # Calculate smear factor
          response_bin_breaks = distinct(data.frame(train_predictions) %>% dplyr::select(c(outer_fold, response_bin_breaks))) # Set up look-up table for response breaks
          correction_factor = distinct(data.frame(train_predictions) %>% dplyr::select(c(outer_fold, correction_factor, response_bin))) # Set up look-up table for correction factors
          
          # Join correction factor to predictions on test set
          final_predictions_subset = left_join(final_predictions_subset, response_bin_breaks, by = 'outer_fold') # Assign response breaks
          final_predictions_subset = final_predictions_subset %>%
            group_by(outer_fold) %>%
            mutate(response_bin = cut(predicted, breaks = response_bin_breaks[[1]], include.lowest = TRUE)) # Assign bin based on breaks from training data
          final_predictions_subset = left_join(data.frame(final_predictions_subset), correction_factor, by = c('outer_fold', 'response_bin')) # Assign correction factor based on outer fold and response bin
          
          # Save uncorrected and corrected predictions
          final_predictions_subset$predicted_bias_adj = exp(final_predictions_subset$predicted) * final_predictions_subset$correction_factor
          final_predictions_subset$predicted = exp(final_predictions_subset$predicted)
        
      }
      
      print('Bias correction applied')
      
      # 2.4 Calculate accuracy metrics - untransformed  ------------------------
      
      # Negative predictions are forced to zero
      final_predictions_subset$predicted[final_predictions_subset$predicted < 0] = 0
      
      # Calculate accuracy metrics
      fit_rmse = round(yardstick::rmse_vec(final_predictions_subset$biomass_density_gm2, final_predictions_subset$predicted), 2)
      fit_mae = round(yardstick::mae_vec(final_predictions_subset$biomass_density_gm2, final_predictions_subset$predicted), 2)
      fit_bias = round(bias_vec(final_predictions_subset$biomass_density_gm2, final_predictions_subset$predicted), 2)
      fit_rsq = round(yardstick::rsq_vec(final_predictions_subset$biomass_density_gm2, final_predictions_subset$predicted), 2)
      fit_mape = round(yardstick::mape_vec(final_predictions_subset$biomass_density_gm2, final_predictions_subset$predicted), 0)
      fit_rmspe = round(rmspe_vec(final_predictions_subset$biomass_density_gm2, final_predictions_subset$predicted), 0)
      
      # Create version of dataset with log-log transformed results
      final_predictions_subset_log = final_predictions_subset
      final_predictions_subset_log$biomass_density_gm2[final_predictions_subset_log$biomass_density_gm2 <= 0] = 1 # Values <= 0 are forced to one to allow log transformation
      final_predictions_subset_log$predicted[final_predictions_subset_log$predicted <= 0] = 1 # Values <= 0 are forced to one to allow log transformation
      final_predictions_subset_log$biomass_density_gm2 = log(final_predictions_subset_log$biomass_density_gm2)
      final_predictions_subset_log$predicted = log(final_predictions_subset_log$predicted)
      
      # Calculate accuracy metrics - log-log
      fit_rmse_log = round(yardstick::rmse_vec(final_predictions_subset_log$biomass_density_gm2, final_predictions_subset_log$predicted), 2)
      fit_mae_log = round(yardstick::mae_vec(final_predictions_subset_log$biomass_density_gm2, final_predictions_subset_log$predicted), 2)
      fit_bias_log = round(bias_vec(final_predictions_subset_log$biomass_density_gm2, final_predictions_subset_log$predicted), 2)
      fit_rsq_log = round(yardstick::rsq_vec(final_predictions_subset_log$biomass_density_gm2, final_predictions_subset_log$predicted), 2)
      fit_mape_log = round(yardstick::mape_vec(final_predictions_subset_log$biomass_density_gm2, final_predictions_subset_log$predicted), 0)
      fit_rmspe_log = round(rmspe_vec(final_predictions_subset_log$biomass_density_gm2, final_predictions_subset_log$predicted), 0)
      
      print('Uncorrected evaluation metrics calculated')
      
      # 2.5 Calculate accuracy metrics - log transformed  ----------------------
      
      if(grepl('log', model_type)){
        
        # Negative predictions are forced to zero
        final_predictions_subset$predicted_bias_adj[final_predictions_subset$predicted_bias_adj < 0] = 0

        # Calculate accuracy metrics
        fit_rmse_bias_adj = round(yardstick::rmse_vec(final_predictions_subset$biomass_density_gm2, final_predictions_subset$predicted_bias_adj), 2)
        fit_mae_bias_adj = round(yardstick::mae_vec(final_predictions_subset$biomass_density_gm2, final_predictions_subset$predicted_bias_adj), 2)
        fit_bias_bias_adj = round(bias_vec(final_predictions_subset$biomass_density_gm2, final_predictions_subset$predicted_bias_adj), 2)
        fit_rsq_bias_adj = round(yardstick::rsq_vec(final_predictions_subset$biomass_density_gm2, final_predictions_subset$predicted_bias_adj), 2)
        fit_mape_bias_adj = round(yardstick::mape_vec(final_predictions_subset$biomass_density_gm2, final_predictions_subset$predicted_bias_adj), 0)
        fit_rmspe_bias_adj = round(rmspe_vec(final_predictions_subset$biomass_density_gm2, final_predictions_subset$predicted_bias_adj), 0)
        
        # Create version of dataset with log-log transformed results
        final_predictions_subset_log_bias_adj = final_predictions_subset
        final_predictions_subset_log_bias_adj$biomass_density_gm2[final_predictions_subset_log$biomass_density_gm2 <= 0] = 1 # Values <= 0 are forced to one to allow log transformation
        final_predictions_subset_log_bias_adj$predicted_bias_adj[final_predictions_subset_log$predicted_bias_adj <= 0] = 1 # Values <= 0 are forced to one to allow log transformation
        final_predictions_subset_log_bias_adj$biomass_density_gm2 = log(final_predictions_subset_log$biomass_density_gm2)
        final_predictions_subset_log_bias_adj$predicted_bias_adj = log(final_predictions_subset_log$predicted_bias_adj)
        
        # Calculate accuracy metrics - log-log
        fit_rmse_log_bias_adj = round(yardstick::rmse_vec(final_predictions_subset_log$biomass_density_gm2, final_predictions_subset_log$predicted_bias_adj), 2)
        fit_mae_log_bias_adj = round(yardstick::mae_vec(final_predictions_subset_log$biomass_density_gm2, final_predictions_subset_log$predicted_bias_adj), 2)
        fit_bias_log_bias_adj = round(bias_vec(final_predictions_subset_log$biomass_density_gm2, final_predictions_subset_log$predicted_bias_adj), 2)
        fit_rsq_log_bias_adj = round(yardstick::rsq_vec(final_predictions_subset_log$biomass_density_gm2, final_predictions_subset_log$predicted_bias_adj), 2)
        fit_mape_log_bias_adj = round(yardstick::mape_vec(final_predictions_subset_log$biomass_density_gm2, final_predictions_subset_log$predicted_bias_adj), 0)
        fit_rmspe_log_bias_adj = round(rmspe_vec(final_predictions_subset_log$biomass_density_gm2, final_predictions_subset_log$predicted_bias_adj), 0)
        
      }
      
      print('Bias corrected evaluation metrics calculated (if appropriate)')
      
      # 2.6 Tidy and combine - untransformed -----------------------------------
      
      mod_original = data.frame(ds_type = ds_type, 
                                mod = model_type,
                                ds_subset = ds_subset,
                                axes = 'original',
                                bias_adj = 'uncorrected',
                                rmse = fit_rmse,
                                mae = fit_mae,
                                bias = fit_bias,
                                rsq = fit_rsq,
                                mape = fit_mape,
                                rmspe = fit_rmspe)
      
      mod_log_log = data.frame(ds_type = ds_type, 
                               mod = model_type,
                               ds_subset = ds_subset,
                               axes = 'log-log',
                               bias_adj = 'uncorrected',
                               rmse = fit_rmse_log,
                               mae = fit_mae_log,
                               bias = fit_bias_log,
                               rsq = fit_rsq_log,
                               mape = fit_mape_log,
                               rmspe = fit_rmspe_log)
      
      mod_uncorrected = rbind(mod_original, mod_log_log)
      
      print('Uncorrected evaluation metrics gathered')
      
      # 2.6 Tidy and combine - log transformed ---------------------------------
      
      if(grepl('log', model_type)){
      
        mod_original_bias_adj = data.frame(ds_type = ds_type, 
                                           mod = model_type,
                                           ds_subset = ds_subset,
                                           axes = 'original',
                                           bias_adj = 'binned_smear',
                                           rmse = fit_rmse_bias_adj,
                                           mae = fit_mae_bias_adj,
                                           bias = fit_bias_bias_adj,
                                           rsq = fit_rsq_bias_adj,
                                           mape = fit_mape_bias_adj,
                                           rmspe = fit_rmspe_bias_adj)
        
        mod_log_log_bias_adj = data.frame(ds_type = ds_type, 
                                          mod = model_type,
                                          ds_subset = ds_subset,
                                          axes = 'log-log',
                                          bias_adj = 'binned_smear',
                                          rmse = fit_rmse_log_bias_adj,
                                          mae = fit_mae_log_bias_adj,
                                          bias = fit_bias_log_bias_adj,
                                          rsq = fit_rsq_log_bias_adj,
                                          mape = fit_mape_log_bias_adj,
                                          rmspe = fit_rmspe_log_bias_adj)
        
        mod_bias_adj = rbind(mod_original_bias_adj, mod_log_log_bias_adj)
        
        mod_df = rbind(mod_uncorrected, mod_bias_adj)
        
      }else{
        
        mod_df = mod_uncorrected
        
      }
      
      print('Bias corrected evaluation metrics gathered (if appropriate)')
      
      mod_evaluation = rbind(mod_evaluation, mod_df)

      print('Final evaluation metrics added to dataset')
      cat('\n')
      
    } # END ZEROS HANDLING LOOP
    
  } # END MODEL TYPE LOOP
  
} # END DATASET TYPE LOOP

# ==============================================================================
# 3. SAVE ======================================================================
# ==============================================================================

write.csv(mod_evaluation, paste0(out_dir, 'continuous_model_accuracy_comparison.csv'), row.names = FALSE)


