################################################################################
################################################################################

# DESCRIPTION:

# This script produces plots of final model accuracy
# - Combines classification + regression
# - Uses all data, including zeros
# - For log transformed models, bias corrected predictions are used

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
library(probably)
library(ggplot2)
library(ggpmisc)
library(yardstick)
library(classInt)
library(viridis)
library(cowplot)
library(grid)
library(MethComp)

source('scripts/000.0_yardstick_bias.R')
source('scripts/000.1_yardstick_rmspe.R')

# 1.2 Parameters ---------------------------------------------------------------

ds_types = c('total', 'woody')
model_binary = 'rf'
model_continuous = 'rf'
binary_threshold_metric = 'j_index' # Metric used to choose best threshold
bin_max = 200 # Choose based on automatically selected bins, but specify here to remain consistent across plots - for 'all' data should be around 200, for '' data should be around

dir_binary = paste0('output/08_cv_results/binary_loocv_10/')
dir_continuous = paste0('output/08_cv_results/continuous_loocv_10/')
dir_out = paste0('output/10_cv_final_model/')

# ==============================================================================
# 2. CREATE FUNCTION ===========================================================
# ==============================================================================

make_plots = function(ds_type){ # START FUNCTION
  
  # ==============================================================================
  # 3. GATHER FILES AND TIDY DATASETS ============================================
  # ==============================================================================
  
  # 3.1 Get files --------------------------------------------------------------
  
  threshold_files = list.files(dir_binary, paste0('classification_thresholds_', ds_type, '_', model_binary), full.names = TRUE)
  predictions_binary_files = list.files(dir_binary, paste0('final_predictions_', ds_type, '_', model_binary), full.names = TRUE)
  predictions_continuous_files = list.files(dir_continuous, paste0('final_predictions_', ds_type, '_', model_continuous), full.names = TRUE)
  train_predictions_files = list.files(dir_continuous, paste0('final_predictions_train_', ds_type, '_', model_continuous), full.names = TRUE)
  
  # Remove log transformed models if not specified
  if(!grepl('log', model_binary)){
    threshold_files = threshold_files[!grepl('_log', threshold_files)]
    predictions_binary_files = predictions_binary_files[!grepl('_log', predictions_binary_files)]
  }
  if(!grepl('log', model_continuous)){
    predictions_continuous_files = predictions_continuous_files[!grepl('_log', predictions_continuous_files)]
    train_predictions_files = train_predictions_files[!grepl('_log', train_predictions_files)]
  }
  
  # 3.2 Aggregate data ---------------------------------------------------------
  
  thresholds_binary = do.call(rbind, lapply(threshold_files, read.csv, header = TRUE))
  predictions_binary = do.call(rbind, lapply(predictions_binary_files, read.csv, header = TRUE))
  predictions_continuous = do.call(rbind, lapply(predictions_continuous_files, read.csv, header = TRUE))
  train_predictions = do.call(rbind, lapply(train_predictions_files, read.csv, header = TRUE))
  
  # ============================================================================
  # 4. ASSIGN CLASSIFICATIONS ==================================================
  # ============================================================================
  
  # 4.1 Order classification variable properly ---------------------------------
  
  predictions_binary$presence_absence = as.factor(predictions_binary$presence_absence)
  predictions_binary$presence_absence = factor(predictions_binary$presence_absence, levels = c('presence', 'absence'))
  
  # 4.2 Choose presence threshold ----------------------------------------------
  
  # Set up F measure metric with beta parameter
  f_meas_beta2 = metric_tweak("f_meas_beta2", f_meas, beta = 2)
  
  # Using predicted probabilities from each outer fold test set, calculate accuracy metrics anew for each potential presence threshold    
  thresholds_binary = predictions_binary %>%
    mutate(presence_absence = as.factor(presence_absence)) %>%
    threshold_perf(truth = presence_absence, 
                   estimate = prob_presence, 
                   thresholds = seq(0, 1, by = 0.05),
                   metrics = metric_set(yardstick::accuracy, yardstick::bal_accuracy, yardstick::kap, yardstick::f_meas, f_meas_beta2, yardstick::j_index, yardstick::mcc, yardstick::sens, yardstick::spec, yardstick::precision))
  names(thresholds_binary) = gsub('[.]', '', names(thresholds_binary))
  
  # Get best threshold from calculated results
  best_threshold = thresholds_binary[thresholds_binary$metric == binary_threshold_metric,] %>%
    group_by(threshold) %>%
    summarise(estimate_avg = mean(estimate)) %>%
    slice_max(estimate_avg, n = 1) %>% # Choose threshold with maximum metric (j_index, sensitivity or specificity)
    slice_min(threshold, n = 1) # If there is a tie, choose the lower threshold
  
  # 4.3 Assign classifications -------------------------------------------------
  
  predictions_binary = predictions_binary %>%
    mutate(
      predicted = make_two_class_pred(
        levels = levels(presence_absence), # A character vector of class levels
        estimate = prob_presence, # A single numeric vector corresponding to the class probabilities of the first level in levels
        threshold = best_threshold$threshold
      )
    )
  
  # ============================================================================
  # 5. APPLY BACK-TRANSFORMATION BIAS CORRECTION IF SPECIFIED ==================
  # ============================================================================
  
  if(grepl('log', model_continuous)){
    
    if(bias_adj_binned){
      
      # Calculate correction factor - binned smear
      train_predictions = train_predictions %>%
        group_by(outer_fold) %>%
        mutate(response_bin_breaks = list(c(head(log(classIntervals(exp(predicted), n=5, style="fisher", warnLargeN = FALSE)$brks), -1), Inf))) %>% # Calculates breaks, replaces last value with Inf to accommodate values in test predictions that exceed those in training predictions
        mutate(response_bin = cut(predicted, breaks = response_bin_breaks[[1]], include.lowest = TRUE)) # Assigns bins based on breaks
      train_predictions = train_predictions %>%
        group_by(outer_fold, response_bin) %>% # Group by outer fold AND response bin
        mutate(correction_factor = sum(exp(biomass_density_gm2  - predicted))/length(biomass_density_gm2)) # Calculate smear factor
      response_bin_breaks = distinct(data.frame(train_predictions) %>% select(c(outer_fold, response_bin_breaks))) # Set up look-up table for response breaks
      correction_factor = distinct(data.frame(train_predictions) %>% select(c(outer_fold, correction_factor, response_bin))) # Set up look-up table for correction factors
      
      # Join correction factor to predictions on test set
      predictions_continuous = left_join(predictions_continuous, response_bin_breaks, by = 'outer_fold') # Assign response breaks
      predictions_continuous = predictions_continuous %>%
        group_by(outer_fold) %>%
        mutate(response_bin = cut(predicted, breaks = response_bin_breaks[[1]], include.lowest = TRUE)) # Assign bin based on breaks from training data
      predictions_continuous = left_join(data.frame(predictions_continuous), correction_factor, by = c('outer_fold', 'response_bin')) # Assign correction factor based on outer fold and response bin
      
    }else{
      
      # Calculate correction factor
      train_predictions = train_predictions %>%
        group_by(outer_fold) %>%
        mutate(correction_factor = lm(exp(biomass_density_gm2) ~ 0 + exp(predicted))$coefficients[[1]]) # Match moments
      correction_factor = distinct(train_predictions %>% select(c(outer_fold, correction_factor))) # Set up look-up table for correction factors
      
      # Join correction factor to predictions on test set
      predictions_continuous = left_join(predictions_continuous, correction_factor, by = 'outer_fold') # Assign correction factor based on outer fold
      
    }
    
    if(bias_adj){
      
      predictions_continuous$predicted = exp(predictions_continuous$predicted) * predictions_continuous$correction_factor
      
    }else{
      
      predictions_continuous$predicted = exp(predictions_continuous$predicted)
      
    }
    
  }
  
  # ============================================================================
  # 6. COMBINE CLASSIFICATION AND REGRESSION RESULTS ===========================
  # ============================================================================
  
  # 6.1 Join -------------------------------------------------------------------
  
  predictions_binary = predictions_binary %>% rename(predicted_binary = predicted)
  predictions = left_join(predictions_continuous, predictions_binary %>% dplyr::select(c(site_code, predicted_binary, presence_absence, prob_presence, prob_absence)))
  
  # 6.2 Handle absences ----------------------------------------------------------
  
  # Negative predictions are forced to zero
  predictions$predicted[predictions$predicted < 0] = 0
  
  # Recode classification
  predictions = predictions %>% mutate(predicted_binary = recode(as.character(predicted_binary), 'presence' = 1, 'absence' = 0))
  
  # Set predicted absences to zero
  predictions$predicted = predictions$predicted * as.numeric(predictions$predicted_binary)
  
  # ============================================================================
  # 7. CALCULATE ACCURACY METRICS ==============================================
  # ============================================================================
  
  # Calculate accuracy metrics
  fit_rmse = round(yardstick::rmse_vec(predictions$biomass_density_gm2, predictions$predicted), 2)
  fit_mae = round(yardstick::mae_vec(predictions$biomass_density_gm2, predictions$predicted), 2)
  fit_rsq = round(yardstick::rsq_vec(predictions$biomass_density_gm2, predictions$predicted), 2)
  fit_bias = round(bias_vec(predictions$biomass_density_gm2, predictions$predicted), 2)
  fit_relative_rmse = round(fit_rmse/mean(predictions$biomass_density_gm2), 2)
  fit_relative_mae = round(fit_mae/mean(predictions$biomass_density_gm2), 2)
  
  # Create log-log transformed dataset
  predictions_log = predictions
  predictions_log$biomass_density_gm2[predictions_log$biomass_density_gm2 < 1] = 1 # Zero and near zero values are forced to one to allow log transformation
  predictions_log$predicted[predictions_log$predicted < 1] = 1 # Zero and near zero values are forced to one to allow log transformation
  predictions_log$biomass_density_gm2 = log(predictions_log$biomass_density_gm2)
  predictions_log$predicted = log(predictions_log$predicted)
  
  # Calculate accuracy metrics - log-log
  fit_rmse_log = round(yardstick::rmse_vec(predictions_log$biomass_density_gm2, predictions_log$predicted), 2)
  fit_mae_log = round(yardstick::mae_vec(predictions_log$biomass_density_gm2, predictions_log$predicted), 2)
  fit_bias_log = round(bias_vec(predictions_log$biomass_density_gm2, predictions_log$predicted), 2)
  fit_rsq_log = round(yardstick::rsq_vec(predictions_log$biomass_density_gm2, predictions_log$predicted), 2)
  fit_relative_rmse_log = round(fit_rmse_log/mean(predictions_log$biomass_density_gm2), 2)
  fit_relative_mae_log = round(fit_mae_log/mean(predictions_log$biomass_density_gm2), 2)
  
  # ============================================================================
  # 8. PLOT ====================================================================
  # ============================================================================
  
  # 8.1 Plot parameters -------------------------------------------------------
  
  # Palettes
  cb_palette_long = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                      "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  
  cb_palette = c("#88CCEE", "#CC6677", "#888888", "#117733")
  
  # Parameters
  biomass_max = max(max(predictions$predicted), max(predictions$biomass_density_gm2))
  textformula = y ~ x
  base_size = 24
  eq_size = base_size/3
  
  # Get title 
  ds_type_title = if(ds_type == 'total'){'Plant'}else if(ds_type == 'woody'){'Woody Plant'}else{stop('ds_type not recognized')}
  
  # 8.2 Plot - normal axes -----------------------------------------------------
  
  # Calculate total least squares regression i.e. orthogonal regression
  tls = Deming(y=predictions$biomass_density_gm2, x=predictions$predicted)
  
  # Plot
  plt_norm = 
    ggplot(predictions, aes(x = predicted,  y = biomass_density_gm2))+
    geom_hex(bins = 15)+
    geom_abline(slope = 1, intercept = 0, lty = 2, size = 1.5)+
    geom_abline(slope = tls['Slope'], intercept = tls['Intercept'], col = 'black', size = 1.5)+
    ylim(c(-1, biomass_max))+
    xlim(c(-1, biomass_max))+
    scale_fill_viridis(option = 'mako', 
                       direction = -1, 
                       limits = c(0, bin_max),
                       oob = scales::squish)+
    labs(x = expression(paste('Biomass Predicted (g ', m^-2, ')')), 
         y = expression(paste('Biomass Observed (g ', m^-2, ')')), 
         col = '', 
         fill = '',
         title = ds_type_title)+
    theme_minimal(base_size = base_size)+
    theme(plot.title = element_text(hjust = 0.5))+
    coord_fixed(ratio = 1)+
    stat_poly_eq(size = eq_size, formula = textformula, geom = 'label', aes(label = paste(after_stat(rr.label), sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1.05, vjust = -1.85, color = 'grey20', fill = 'white', label.size = NA, alpha = 0, label.padding = unit(0.01, "lines"))+
    geom_text(size = eq_size, color = 'grey20', data = predictions[1,], aes(label = paste0('MAE = ', fit_mae)), x = Inf, y = -Inf, hjust = 1.05, vjust = -2.5)+
    geom_text(size = eq_size, color = 'grey20', data = predictions[1,], aes(label = paste0('RMSE = ', fit_rmse)), x = Inf, y = -Inf, hjust = 1.05, vjust = -1)+
    geom_text(size = eq_size, color = 'grey20', data = predictions[1,], aes(label = paste0('rMAE = ', fit_relative_mae)), x = -Inf, y = Inf, hjust = -0.05, vjust = 3.5)+
    geom_text(size = eq_size, color = 'grey20', data = predictions[1,], aes(label = paste0('rRMSE = ', fit_relative_rmse)), x = -Inf, y = Inf, hjust = -0.05, vjust = 2)
  
  # 8.3 Plot - log-log axes ----------------------------------------------------
  
  # Calculate total least squares regression i.e. orthogonal regression
  tls_log = Deming(y=predictions_log$biomass_density_gm2, x=predictions_log$predicted) # TLS i.e. Orthogonal Regression
  
  # Set plot parameters
  biomass_min_log = min(min(predictions_log$predicted), min(predictions_log$biomass_density_gm2))
  biomass_max_log = max(max(predictions_log$predicted), max(predictions_log$biomass_density_gm2))
  
  # Plot
  plt_log_log = 
    ggplot(predictions_log, aes(x = predicted,  y = biomass_density_gm2))+
    ylim(c(biomass_min_log-0.1, biomass_max_log))+
    xlim(c(biomass_min_log-0.1, biomass_max_log))+
    geom_hex(bins = 15)+
    geom_abline(slope = 1, intercept = 0, lty = 2, size = 1.5)+
    geom_abline(slope = tls_log['Slope'], intercept = tls_log['Intercept'], col = 'black', size = 1.5)+
    scale_fill_viridis(option = 'mako', 
                       direction = -1, 
                       limits = c(0, bin_max),
                       oob = scales::squish)+
    labs(x = expression(paste('Biomass Predicted (log(g ', m^-2, '))')), 
         y = expression(paste('Biomass Observed (log(g ', m^-2, '))')), 
         col = '', 
         fill = '',
         title = ds_type_title)+
    theme_minimal(base_size = base_size)+
    theme(plot.title = element_text(hjust = 0.5))+
    coord_fixed(ratio = 1)+
    stat_poly_eq(size = eq_size, formula = textformula, geom = 'label', aes(label = paste(after_stat(rr.label), sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1.05, vjust = -1.85, color = 'grey20', fill = 'white', label.size = NA, alpha = 0, label.padding = unit(0.01, "lines"))+
    geom_text(size = eq_size, color = 'grey20', data = predictions_log[1,], aes(label = paste0('MAE = ', fit_mae_log)), x = Inf, y = -Inf, hjust = 1.05, vjust = -2.5)+
    geom_text(size = eq_size, color = 'grey20', data = predictions_log[1,], aes(label = paste0('RMSE = ', fit_rmse_log)), x = Inf, y = -Inf, hjust = 1.05, vjust = -1)+
    geom_text(size = eq_size, color = 'grey20', data = predictions_log[1,], aes(label = paste0('rMAE = ', fit_relative_mae_log)), x = -Inf, y = Inf, hjust = -0.05, vjust = 3.5)+
    geom_text(size = eq_size, color = 'grey20', data = predictions_log[1,], aes(label = paste0('rRMSE = ', fit_relative_rmse_log)), x = -Inf, y = Inf, hjust = -0.05, vjust = 2)
  
  # 8.4 Create list of plots ---------------------------------------------------
  
  plt_list_inner = list()
  plt_list_inner[[paste0(ds_type, '_norm')]] = plt_norm
  plt_list_inner[[paste0(ds_type, '_log_log')]] = plt_log_log
  
  return(plt_list_inner)
  
} # END FUNCTION

# ==============================================================================
# 9. APPLY FUNCTION ============================================================
# ==============================================================================

plts_all = lapply(ds_types, make_plots)
plt_list = unlist(plts_all, recursive = FALSE)

# ==============================================================================
# 10. COMBINE PLOTS =============================================================
# ==============================================================================

# 10.1 Create supplementary items -----------------------------------------------

base_size = 24

# Legend
legend = cowplot::get_legend(
  plt_list[['total_norm']] +
    theme(legend.position = "right",
          legend.margin = margin(l = 30, b = 100),
          legend.key.height = unit(1.25, "cm"),
          legend.key.width  = unit(0.75, "cm"),
          legend.text = element_text(size = 24))
)

# X-axis text - normal
xtext_norm = ggdraw() +
  draw_label(expression(paste('Biomass Predicted (g ', m^-2, ')')),
             hjust = 0.5,
             vjust = 0,
             size = base_size)

# X-axis text - log-log
xtext_log_log = ggdraw() +
  draw_label(expression(paste('Biomass Predicted (log(g ', m^-2, '))')),
             hjust = 0.5,
             vjust = 0,
             size = base_size)

# Dataset type titles
title_total = plt_list[['total_norm']]
title_total$layers = NULL
title_woody = plt_list[['woody_norm']]
title_woody$layers = NULL

titles = cowplot::plot_grid(title_total + theme_void() + theme(plot.title = element_text(size = base_size + 10)), 
                            NULL, 
                            title_woody + theme_void() + theme(plot.title = element_text(size = base_size + 10)), 
                            nrow = 1,
                            align = 'v',
                            rel_widths = c(1, -0.4, 1))

# 10.2 Create top row - regular axes --------------------------------------------

# Plots
# Blank space hack: https://github.com/wilkelab/cowplot/issues/31
top_row = cowplot::plot_grid(plt_list[['total_norm']] +
                               theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
                                     panel.spacing = unit(c(0, 0, 0, 0), "cm"),
                                     panel.border=element_blank(),
                                     legend.position = "none",
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.x = element_blank(),
                                     plot.title = element_blank()), 
                             NULL, 
                             plt_list[['woody_norm']] +
                               theme(axis.text.y = element_blank(),
                                     axis.ticks.y = element_blank(),
                                     axis.title.y = element_blank(),
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.x = element_blank(),
                                     axis.ticks.length = unit(0, "cm"),
                                     axis.ticks.margin = unit(0, "cm"),
                                     plot.title = element_blank(),
                                     plot.margin = unit(c(0, 0, 0, 0), "cm"),
                                     panel.spacing = unit(c(0, 0, 0, 0), "cm"),
                                     panel.border=element_blank(),
                                     legend.position = "none"), 
                             nrow = 1,
                             align = 'v',
                             labels = c('a)', 'b)'),
                             label_size = base_size + 10,
                             hjust = c(-7, 2.5),
                             vjust = -0.1,
                             rel_widths = c(1, -0.4, 1))

# 10.3 Create bottom row - log-log axes --------------------------------------------

# Plots
# Blank space hack: https://github.com/wilkelab/cowplot/issues/31
bottom_row = cowplot::plot_grid(plt_list[['total_log_log']] +
                                  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
                                        panel.spacing = unit(c(0, 0, 0, 0), "cm"),
                                        panel.border=element_blank(),
                                        legend.position = "none",
                                        axis.text.x = element_blank(),
                                        axis.ticks.x = element_blank(),
                                        axis.title.x = element_blank(),
                                        plot.title = element_blank()), 
                                NULL, 
                                plt_list[['woody_log_log']] +
                                  theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.text.x = element_blank(),
                                        axis.ticks.x = element_blank(),
                                        axis.title.x = element_blank(),
                                        axis.ticks.length = unit(0, "cm"),
                                        axis.ticks.margin = unit(0, "cm"),
                                        plot.title = element_blank(),
                                        plot.margin = unit(c(0, 0, 0, 0), "cm"),
                                        panel.spacing = unit(c(0, 0, 0, 0), "cm"),
                                        panel.border=element_blank(),
                                        legend.position = "none"), 
                                nrow = 1,
                                align = 'v',
                                labels = c('c)', 'd)'),
                                label_size = base_size + 10,
                                hjust = c(-7, 2.5),
                                vjust = -0.1,
                                rel_widths = c(1, -0.4, 1))

# 10.4 Put together --------------------------------------------

# Add text
plt_xtext = cowplot::plot_grid(titles, top_row, xtext_norm, bottom_row, xtext_log_log, rel_heights = c(0.2, 1, 0.25, 1, 0.17), nrow = 5)

# Add legend
plt_final = cowplot::plot_grid(plt_xtext, NULL, legend, rel_widths = c(3, -0.4, 0.4), nrow = 1)

# Save
ggsave(
  paste0(dir_out, 'cv_final_model_accuracy_all_data.png'),
  plt_final,
  width = 40,
  height = 30,
  units = 'cm',
  bg = 'white',
  dpi = 600
)