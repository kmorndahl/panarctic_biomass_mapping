################################################################################
################################################################################

# DESCRIPTION:

# This script summarizes regression performance across model types

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

# NOTES:
# - Final accuracy figures are created using:
#   - All data (including absences)
#   - Bias corrected predictions for log transformed models

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

version = 'v20240514'
ds_types = c('total', 'woody')
model_list = c('rf')
metric_type = 'rmse' # Metric used to determine best hyperparameters
include_obs = FALSE # Include observations on plots (TRUE) or just plot line fits (FALSE)

in_dir = paste0('output/08_cv_results/continuous_loocv_10/')
out_dir_tuning = paste0('output/09_cv_model_comparison/continuous_loocv_10/tuning_plots/')
out_dir_accuracy = paste0('output/09_cv_model_comparison/continuous_loocv_10/accuracy_plots/')

# ==============================================================================
# 2. LOOP DATASETS AND MODELS AND EXPORT FIGURES ===============================
# ==============================================================================

for(ds_type in ds_types){ # START DATASET TYPE LOOP
  
  print(paste0('The dataset type is: ', ds_type))
  cat('\n')
  
  for(model_type in model_list){ # START MODEL TYPE LOOP
    
    print(paste0('The model type is: ', model_type))
    cat('\n')
    
    # 2.1 Gather and combine files ---------------------------------------------
    
    accuracy_metrics_files = list.files(paste0(in_dir, '/accuracy_metrics/'), paste0('accuracy_metrics_', ds_type, '_', model_type), full.names = TRUE)
    train_predictions_files = list.files(paste0(in_dir, '/predictions_train/'), paste0('final_predictions_train_', ds_type, '_', model_type), full.names = TRUE)
    final_predictions_files = list.files(paste0(in_dir, '/predictions/'), paste0('final_predictions_', ds_type, '_', model_type), full.names = TRUE)

    # Remove log transformed models if not specified
    if(!grepl('log', model_type)){
      accuracy_metrics_files = accuracy_metrics_files[!grepl('_log', accuracy_metrics_files)]
      train_predictions_files = train_predictions_files[!grepl('_log', train_predictions_files)]
      final_predictions_files = final_predictions_files[!grepl('_log', final_predictions_files)]
    }
    
    accuracy_metrics = do.call(rbind, lapply(accuracy_metrics_files, read.csv, header = TRUE))
    train_predictions = do.call(rbind, lapply(train_predictions_files, read.csv, header = TRUE))
    final_predictions = do.call(rbind, lapply(final_predictions_files, read.csv, header = TRUE))

    print('Data read in')
      
    # 2.2 Gather and tidy hyperparameter tuning data ---------------------------
    
    # Assign factors
    accuracy_metrics$outer_fold = as.factor(accuracy_metrics$outer_fold)
    accuracy_metrics$inner_fold = as.factor(accuracy_metrics$id)
    accuracy_metrics = accuracy_metrics %>% dplyr::select(-id)
    
    if(model_type != 'linearBaseline'){ # linearBaseline does not have parameters to tune
    
      # Get average metrics for each outer fold, workflow (i.e. # of features) and configuration (i.e. set of parameters)
      # Averages across inner folds
      metrics_avg = accuracy_metrics %>%
        dplyr::select(-c(preproc, model, inner_fold, .estimator, log)) %>%
        group_by(outer_fold, wflow_id, .config, .metric) %>%
        summarise(
          across(.estimate, mean),
          across(-.estimate, first)
        )
    
      # Get best parameter combo per metric type, outer fold and workflow (i.e. # of features)
    
      chosen_params_min = metrics_avg[metrics_avg$.metric %in% c('mae', 'mape', 'rmse', 'rmspe'),] %>%
        group_by(outer_fold, wflow_id, .metric) %>%
        slice_min(.estimate, n=1)
      
      chosen_params_abs_min = metrics_avg[metrics_avg$.metric == 'bias',] %>%
        group_by(outer_fold, wflow_id, .metric) %>%
        slice_min(abs(.estimate), n=1)
      
      chosen_params_max = metrics_avg[metrics_avg$.metric == 'rsq',] %>%
        group_by(outer_fold, wflow_id, .metric) %>%
        slice_max(.estimate, n=1)
      
      chosen_params = rbind(chosen_params_min, chosen_params_abs_min, chosen_params_max)
      
      # Get parameter names
      params = names(data.frame(chosen_params) %>% dplyr::select(-c('outer_fold', 'wflow_id', '.config', '.metric', '.estimate', 'n_feat', 'feat_thresh')))
    
    }
    
    # 2.3 Loop through hyperparameter sets and plot-----------------------------
    
    if(model_type != 'linearBaseline'){ # linearBaseline does not have parameters to tune
    
      for(param in params){
    
        metrics_avg$n_feat = as.factor(metrics_avg$n_feat)
        
        if(include_obs){
          plt = ggplot(metrics_avg[metrics_avg$.metric == metric_type,], aes_string(x = param, y = ".estimate", col = "n_feat"))+
            geom_line(aes(group = interaction(as.factor(n_feat), outer_fold)), alpha = 0.01)+
            geom_point(size = 2, alpha = 0.01)+
            geom_smooth(aes(group = as.factor(n_feat)), se = TRUE, size = 1.5, method = "loess")+
            labs(y = 'RMSE', col = '# of \npredictors')+
            theme_minimal(base_size = 30)
        }else{
          plt = ggplot(metrics_avg[metrics_avg$.metric == metric_type,], aes_string(x = param, y = ".estimate", col = "n_feat"))+
            geom_smooth(aes(group = as.factor(n_feat)), se = TRUE, size = 1.5, method = "loess")+
            labs(y = 'RMSE', col = '# of \npredictors')+
            theme_minimal(base_size = 30)
        }
        
        ggsave(
          paste0(out_dir_tuning, 'parameter_tuning_', ds_type, '_', model_type, '_', param, '.png'),
          plt,
          width = 40,
          height = 30,
          units = 'cm',
          bg = 'white',
          dpi = 150
        )
    
      }
    
    }
    
    # 2.4 Apply bias correction to log transformed models  -------------------
    
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
      final_predictions = left_join(final_predictions, response_bin_breaks, by = 'outer_fold') # Assign response breaks
      final_predictions = final_predictions %>%
        group_by(outer_fold) %>%
        mutate(response_bin = cut(predicted, breaks = response_bin_breaks[[1]], include.lowest = TRUE)) # Assign bin based on breaks from training data
      final_predictions = left_join(data.frame(final_predictions), correction_factor, by = c('outer_fold', 'response_bin')) # Assign correction factor based on outer fold and response bin
  
      # Accuracy plots are created with bias corrected predictions
      final_predictions$predicted = exp(final_predictions$predicted) * final_predictions$correction_factor

    }
    
    # 2.5 Calculate accuracy metrics -------------------------------------------
    
    # Negative predictions are forced to zero
    final_predictions$predicted[final_predictions$predicted < 0] = 0
    
    # Calculate accuracy metrics
    fit_rmse = round(yardstick::rmse_vec(final_predictions$biomass_density_gm2, final_predictions$predicted), 2)
    fit_mae = round(yardstick::mae_vec(final_predictions$biomass_density_gm2, final_predictions$predicted), 2)
    fit_bias = round(bias_vec(final_predictions$biomass_density_gm2, final_predictions$predicted), 2)
    fit_rsq = round(yardstick::rsq_vec(final_predictions$biomass_density_gm2, final_predictions$predicted), 2)
    fit_mape = round(yardstick::mape_vec(final_predictions$biomass_density_gm2, final_predictions$predicted), 0)
    fit_rmspe = round(rmspe_vec(final_predictions$biomass_density_gm2, final_predictions$predicted), 0)
    
    # Create version of dataset with log-log transformed results
    final_predictions_log = final_predictions
    final_predictions_log$biomass_density_gm2[final_predictions_log$biomass_density_gm2 <= 0] = 1 # Values <= 0 are forced to one to allow log transformation
    final_predictions_log$predicted[final_predictions_log$predicted <= 0] = 1 # Values <= 0 are forced to one to allow log transformation
    final_predictions_log$biomass_density_gm2 = log(final_predictions_log$biomass_density_gm2)
    final_predictions_log$predicted = log(final_predictions_log$predicted)
    
    # Calculate accuracy metrics - log-log
    fit_rmse_log = round(yardstick::rmse_vec(final_predictions_log$biomass_density_gm2, final_predictions_log$predicted), 2)
    fit_mae_log = round(yardstick::mae_vec(final_predictions_log$biomass_density_gm2, final_predictions_log$predicted), 2)
    fit_bias_log = round(bias_vec(final_predictions_log$biomass_density_gm2, final_predictions_log$predicted), 2)
    fit_rsq_log = round(yardstick::rsq_vec(final_predictions_log$biomass_density_gm2, final_predictions_log$predicted), 2)
    fit_mape_log = round(yardstick::mape_vec(final_predictions_log$biomass_density_gm2, final_predictions_log$predicted), 0)
    fit_rmspe_log = round(rmspe_vec(final_predictions_log$biomass_density_gm2, final_predictions_log$predicted), 0)
    
    # 2.6 Set plot parameters --------------------------------------------------
    
    cb_palette_long = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                        "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
    
    cb_palette = c("#88CCEE", "#CC6677", "#888888", "#117733")
    
    biomass_max = max(max(final_predictions$predicted), max(final_predictions$biomass_density_gm2))
    textformula = y ~ x
    base_size = 30
    eq_size = 30/4
    
    # 2.7 Plot - normal axes ---------------------------------------------------
    
    plt = 
      ggplot(final_predictions, aes(x = predicted,  y = biomass_density_gm2))+
      geom_point(aes(col = zone_3, fill = zone_3), shape = 21, col = 'black', size = 4, alpha = 0.75)+
      geom_smooth(method = 'lm', col = 'black', size = 0.5)+
      geom_abline(slope = 1, intercept = 0, lty = 2, size = 0.5)+
      ylim(c(0, biomass_max))+
      xlim(c(0, biomass_max))+
      scale_fill_manual(values = cb_palette)+
      scale_colour_manual(values = cb_palette)+
      labs(x = expression(paste('Total Biomass (g  ', m^-2, ') Predicted')), y = expression(paste('Total Biomass (g  ', m^-2, ') Observed')), fill = 'Locale')+
      theme_minimal(base_size = base_size)+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions[1,], aes(label = paste0('RMSE = ', fit_rmse)), x = Inf, y = -Inf, hjust = 1.05, vjust = -9.5)+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions[1,], aes(label = paste0('MAE = ', fit_mae)), x = Inf, y = -Inf, hjust = 1.05, vjust = -8)+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions[1,], aes(label = paste0('RMSPE = ', fit_rmspe, '%')), x = Inf, y = -Inf, hjust = 1.05, vjust = -6.5)+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions[1,], aes(label = paste0('MAPE = ', fit_mape, '%')), x = Inf, y = -Inf, hjust = 1.05, vjust = -5)+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions[1,], aes(label = paste0('Bias = ', fit_bias)), x = Inf, y = -Inf, hjust = 1.05, vjust = -3.5)+
      stat_poly_eq(size = eq_size, formula = textformula, geom = 'label', aes(label = paste(after_stat(rr.label), sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1.1, vjust = -1, color = 'grey20', fill = 'white', label.size = NA, alpha = 0, label.padding = unit(0.01, "lines"))+
      stat_poly_eq(size = eq_size, formula = textformula, geom = 'label', aes(label = paste(after_stat(eq.label), sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1.05, vjust = -0.3, color = 'grey20', fill = 'white', label.size = NA, alpha = 0, label.padding = unit(0.01, "lines"))
    plt
    
    ggsave(
      paste0(out_dir_accuracy, 'final_assessment_', ds_type, '_', model_type, '.png'),
      plt,
      width = 40,
      height = 30,
      units = 'cm',
      bg = 'white',
      dpi = 150
    )
    
    # 2.8 Plot - log-log axes --------------------------------------------------
    
    biomass_max_log = max(max(final_predictions_log$predicted), max(final_predictions_log$biomass_density_gm2))
    
    plt_log = 
      ggplot(final_predictions_log, aes(x = predicted,  y = biomass_density_gm2))+
      geom_point(aes(col = zone_3, fill = zone_3), shape = 21, col = 'black', size = 4, alpha = 0.75)+
      geom_smooth(method = 'lm', col = 'black', size = 0.5)+
      geom_abline(slope = 1, intercept = 0, lty = 2, size = 0.5)+
      ylim(c(0, biomass_max_log))+
      xlim(c(0, biomass_max_log))+
      scale_fill_manual(values = cb_palette)+
      scale_colour_manual(values = cb_palette)+
      labs(x = expression(paste('Total Biomass (g  ', m^-2, ') Predicted')), y = expression(paste('Total Biomass (g  ', m^-2, ') Observed')), fill = 'Locale')+
      theme_minimal(base_size = base_size)+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions_log[1,], aes(label = paste0('RMSE = ', fit_rmse_log)), x = Inf, y = -Inf, hjust = 1.05, vjust = -9.5)+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions_log[1,], aes(label = paste0('MAE = ', fit_mae_log)), x = Inf, y = -Inf, hjust = 1.05, vjust = -8)+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions_log[1,], aes(label = paste0('RMSPE = ', fit_rmspe_log, '%')), x = Inf, y = -Inf, hjust = 1.05, vjust = -6.5)+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions_log[1,], aes(label = paste0('MAPE = ', fit_mape_log, '%')), x = Inf, y = -Inf, hjust = 1.05, vjust = -5)+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions_log[1,], aes(label = paste0('Bias = ', fit_bias_log)), x = Inf, y = -Inf, hjust = 1.05, vjust = -3.5)+
      stat_poly_eq(size = eq_size, formula = textformula, geom = 'label', aes(label = paste(after_stat(rr.label), sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1.1, vjust = -1, color = 'grey20', fill = 'white', label.size = NA, alpha = 0, label.padding = unit(0.01, "lines"))+
      stat_poly_eq(size = eq_size, formula = textformula, geom = 'label', aes(label = paste(after_stat(eq.label), sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1.05, vjust = -0.3, color = 'grey20', fill = 'white', label.size = NA, alpha = 0, label.padding = unit(0.01, "lines"))
    plt_log
    
    ggsave(
      paste0(out_dir_accuracy, 'final_assessment_log_', ds_type, '_', model_type, '.png'),
      plt_log,
      width = 40,
      height = 30,
      units = 'cm',
      bg = 'white',
      dpi = 150
    )

  } # END MODEL TYPE LOOP
  
} # END DATASET TYPE LOOP

