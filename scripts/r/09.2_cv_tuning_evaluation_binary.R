################################################################################
################################################################################

# DESCRIPTION:

# This script produces plots of classification accuracy

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
model_list = c('rf')
metric_type = 'roc_auc' # Metric used to determine best hyperparameters
binary_threshold_metric = 'j_index' # Metric used to choose best threshold
threshold_deterimination_method = 'calculate' # Choose 'average' or 'calculate'
include_obs = FALSE # Include observations on plots (TRUE) or just plot line fits (FALSE)

in_dir = paste0('output/08_cv_results/binary_loocv_10/')
out_dir_tuning = paste0('output/09_cv_model_comparison/binary_loocv_10/tuning_plots/')
out_dir_accuracy = paste0('output/09_cv_model_comparison/binary_loocv_10/accuracy_plots/')

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
    
    # Get lists of all files
    accuracy_metrics_files = list.files(paste0(in_dir, '/accuracy_metrics/'), paste0('accuracy_metrics_', ds_type, '_', model_type), full.names = TRUE)
    threshold_files = list.files(paste0(in_dir, '/thresholds/'), paste0('classification_thresholds_', ds_type, '_', model_type), full.names = TRUE)
    final_predictions_files = list.files(paste0(in_dir, '/predictions/'), paste0('final_predictions_', ds_type, '_', model_type), full.names = TRUE)

    # Read in as dataframes and combine
    accuracy_metrics = do.call(rbind, lapply(accuracy_metrics_files, read.csv, header = TRUE))
    thresholds = do.call(rbind, lapply(threshold_files, read.csv, header = TRUE))
    final_predictions = do.call(rbind, lapply(final_predictions_files, read.csv, header = TRUE))

    # 2.2 Gather and tidy hyperparameter tuning data ---------------------------
    
    # Assign factors
    accuracy_metrics$outer_fold = as.factor(accuracy_metrics$outer_fold)
    accuracy_metrics$inner_fold = as.factor(accuracy_metrics$id)
    accuracy_metrics = accuracy_metrics %>% dplyr::select(-id)
    
    # Get average metrics for each outer fold, workflow (i.e. # of features) and configuration (i.e. set of hyperparameters)
    # Averages across inner folds
    metrics_avg = accuracy_metrics %>%
      dplyr::select(-c(preproc, model, inner_fold, .estimator)) %>%
      group_by(outer_fold, wflow_id, .config, .metric) %>%
      summarise(
        across(.estimate, mean),
        across(-.estimate, first)
      )
    
    # Get best parameter combo per metric type, outer fold and workflow (i.e. # of features)
    chosen_params = metrics_avg %>%
      group_by(outer_fold, wflow_id, .metric) %>%
      slice_max(.estimate, n=1)
    
    # Get parameter names
    params = names(data.frame(chosen_params) %>% dplyr::select(-c('outer_fold', 'wflow_id', '.config', '.metric', '.estimate', 'n_feat', 'feat_thresh')))
    
    # 2.3 Loop through hyperparameter sets and plot-----------------------------
    
    for(param in params){
      
      metrics_avg$n_feat = as.factor(metrics_avg$n_feat)
      
      if(include_obs){
        plt = ggplot(metrics_avg[metrics_avg$.metric == metric_type,], aes_string(x = param, y = ".estimate", col = "n_feat"))+
          geom_line(aes(group = interaction(as.factor(n_feat), outer_fold)), alpha = 0.01)+
          geom_point(size = 2, alpha = 0.01)+
          geom_smooth(aes(group = as.factor(n_feat)), se = TRUE, size = 1.5, method = "loess")+
          labs(y = 'ROC AUC', col = '# of \npredictors')+
          theme_minimal(base_size = 30)
      }else{
        plt = ggplot(metrics_avg[metrics_avg$.metric == metric_type,], aes_string(x = param, y = ".estimate", col = "n_feat"))+
          geom_smooth(aes(group = as.factor(n_feat)), se = TRUE, size = 1.5, method = "loess")+
          labs(y = 'ROC AUC', col = '# of \npredictors')+
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
    
    # 2.4 Choose presence threshold --------------------------------------------
    
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
    
    # 2.5 Calculate accuracy metrics -------------------------------------------
    
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
    
    # 2.6 Plot model accuracy --------------------------------------------------
    
    # Set plot parameters
    textformula = y ~ x
    base_size = 60
    eq_size = base_size/3
    
    # Calculate presence threshold
    presence_threshold = best_threshold$threshold
    
    # Scatter plot
    plt =
      ggplot()+
      geom_rect(aes(xmin = -Inf, xmax = presence_threshold, ymin = -Inf, ymax = Inf), fill = "#00BFC4", color = NA, alpha = 0.15)+
      geom_rect(aes(xmin = presence_threshold, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "#F8766D", color = NA, alpha = 0.15)+
      geom_point(data = final_predictions, aes(x = prob_presence, y = biomass_density_gm2, col = presence_absence), size = 12, alpha = 0.7)+
      geom_vline(xintercept = presence_threshold, size = 1.5)+
      scale_x_continuous(breaks = seq(0, 1, 0.1))+
      labs(x = 'Predicted Probability of Presence', y = expression(paste('Biomass Observed (g ', m^-2, ')')), col = 'Observed Classification')+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions[1,], aes(label = paste0('Accuracy = ', fit_accuracy)), x = -Inf, y = Inf, hjust = -0.1, vjust = 2)+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions[1,], aes(label = paste0('F-Score = ', fit_f_meas)), x = -Inf, y = Inf, hjust = -0.1, vjust = 3.5)+
      geom_text(size = eq_size, color = 'grey20', data = final_predictions[1,], aes(label = paste0('J-Index = ', fit_j_index)), x = -Inf, y = Inf, hjust = -0.1, vjust = 5)+
      theme_minimal(base_size = base_size)+
      theme(legend.position="top",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.text = element_text(size = 30),
            legend.title = element_text(size = 40))
    plt
    
    # Save
    ggsave(
      paste0(out_dir_accuracy, 'assessment_', ds_type, '_', model_type, '.png'),
      plt,
      width = 40,
      height = 40,
      units = 'cm',
      bg = 'white',
      dpi = 150
    )
    
  } # END MODEL TYPE LOOP
  
} # END DATASET TYPE LOOP


