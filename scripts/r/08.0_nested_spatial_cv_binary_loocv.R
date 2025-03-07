################################################################################
################################################################################

# DESCRIPTION:

# This script runs probability model cross-validation (total and woody):
# - data grouped by spatial blocks
# - data stratified by presence/absence to prevent class imbalance
# - nested cross-validation with 10 inner folds, leave-one-out outer folds
# - random forest model, tuned mtry, max.depth
# - MRMR feature selection
# - tuned classification threshold

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

# REFERENCES:

# - Cross-validation best practices: https://hastie.su.domains/MOOC-Slides/cv_boot.pdf
# - Nested CV example: https://www.tidymodels.org/learn/work/nested-resampling/
# - Recommended pre-processing: https://www.tmwr.org/pre-proc-table.html
# - All recipe step functions: https://recipes.tidymodels.org/reference/index.html
# - Hurdle model example: https://geoffruddock.com/building-a-hurdle-regression-estimator-in-scikit-learn/
# - Workflow sets: https://www.tmwr.org/workflow-sets.html
# - Classification metrics: https://machinelearningmastery.com/threshold-moving-for-imbalanced-classification/

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)

library(rsample)
library(parsnip)
library(recipes)
library(workflows)
library(tune)
library(yardstick)
library(broom)
library(dials)
library(workflowsets)
library(probably)

library(sf)

library(spatialsample)
library(praznik)
library(colino)

library(ranger)
library(xgboost) 
library(kernlab)
library(glmnet)
library(lightgbm)
library(bonsai)

# 1.2 Parameters ---------------------------------------------------------------

# Bash
params = commandArgs(trailingOnly = TRUE)
ds_type = params[1]
outer_fold_num = params[2]

# R
in_version = 'v20240508'
out_version = 'v20240514'
model_list = c('lgbm', 'rf', 'svmLinear', 'svmPoly', 'xgb')
models_normalize = c('svmLinear', 'svmPoly')
presence_thresh = 10 # Presence threshold
metrics = metric_set(yardstick::accuracy, yardstick::kap, yardstick::f_meas, yardstick::j_index, yardstick::roc_auc, yardstick::pr_auc, yardstick::brier_class, yardstick::sens, yardstick::spec, yardstick::mn_log_loss)
grid_n = 20 # Hyperparameter tuning grid size
n_outer_folds = NULL # CV outer fold number - NULL does leave one group out CV
n_inner_folds = 10 # CV inner fold numger
feat_thres_list = c(0.50, 0.75, 0.90, 0.95) # Feature threshold options
selection_metric = 'roc_auc' # Metric for choosing best parameters - choose 'rmse' or 'mae'
group_var = 'spatial_block_id' # CV grouping variable
strata_var = 'absence_in_block' # CV stratification variable

# Directories
in_dir = 'output/06_model_ready_data/'
out_dir = paste0('output/08_cv_results/binary_loocv_10/')

# 1.3 Read in data -------------------------------------------------------------

biomass_df = read.csv(paste0(in_dir, 'ds_', ds_type, '_', in_version, '.csv'))

# 1.4 Report -------------------------------------------------------------------

cat('\n')
print(paste0('Outer fold: ', outer_fold_num))
print(paste0('In version: ', in_dir))
print(paste0('Out version: ', out_dir))
print(paste0('Data set type: ', ds_type))
print(paste0('Models: ', paste(model_list, collapse = ', ')))
print(paste0('Presence threshold: ', presence_thresh))
print(paste0('Hyperparameter grid size: ', grid_n))
print(paste0('Number of outer folds: ', n_outer_folds))
print(paste0('Number of inner folds: ', n_inner_folds))
print(paste0('Feature threshold options: ', paste(feat_thres_list, collapse = ', ')))
print(paste0('Selection metric: ', selection_metric))
print(paste0('Grouping variable: ', group_var))
print(paste0('Stratification variable: ', strata_var))
cat('\n')

# ==============================================================================
# 2. PREPARE DATA ==============================================================
# ==============================================================================

# 2.1 Set up response variable -------------------------------------------------
# There is no common convention on which factor level should automatically be 
# considered the "event" or "positive" result when computing binary classification metrics 
# In yardstick, the default is to use the first level
# Therefore we order our factors so that 'presence' (the event or positive result) is first

biomass_df$presence_absence = 'absence'
biomass_df$presence_absence[biomass_df$biomass_density_gm2 > presence_thresh] = 'presence'
biomass_df$presence_absence = factor(biomass_df$presence_absence, levels = c('presence', 'absence'))

# 2.2 Prepare spatial data -----------------------------------------------------

# Get lat/long
latlong = biomass_df %>% select(any_of(c("longitude", "latitude")))

# Convert data to sf
biomass_df = sf::st_as_sf(biomass_df, coords = c("longitude", "latitude"), crs = 4326)
biomass_df = cbind(biomass_df, latlong)

# Transform biomass data
biomass_df = st_transform(biomass_df, 3571)

# 2.3 Tidy predictors -----------------------------------------------------

# Specify factors
biomass_df = biomass_df %>% mutate_at(c('zone_3', 'ecoregion', 'land_cover', 'world_terrestrial_ecosystems'), as.factor)

# Initialize full list of predictors
predictors = grep("permafrost|spectral_|texture_|topo_|trees_|zone_3|ecoregion|world_terrestrial_ecosystems|latitude|longitude", names(biomass_df), value = TRUE)

# 2.4 Set up spatial blocks ----------------------------------------------------

set.seed(1908)
block_folds = spatialsample::spatial_block_cv(biomass_df, v = NULL, cellsize = 100*1000) # Cellsize units should be in m: 100 km times 1000 m in a km

print('Spatial block folds created:')
print(block_folds)
cat('\n')

# 2.5 Get spatial block IDs ----------------------------------------------------

# Fetch spatial block IDs
spatial_block_df = data.frame()

for(b_row in 1:nrow(block_folds)){
  
  fold = block_folds[b_row,]
  test = assessment(fold$splits[[1]])
  block_id = gsub('Fold', '', fold$id)
  
  spatial_block_df = rbind(spatial_block_df, cbind(test, data.frame('spatial_block_id' = block_id)))
  
}

# Join spatial block ID to main data set
spatial_block_df = spatial_block_df %>% dplyr::select(site_code, spatial_block_id)
spatial_block_df$spatial_block_id = as.factor(spatial_block_df$spatial_block_id)
spatial_block_df = st_drop_geometry(spatial_block_df)
biomass_df = st_drop_geometry(biomass_df)
biomass_df = left_join(biomass_df, spatial_block_df)

# 2.6 Set up stratification variable -------------------------------------------

# Assign 'nonveg' to a spatial block if any observation within that block has a nonveg site
biomass_df = biomass_df %>%
  group_by(spatial_block_id) %>%
  mutate(absence_in_block = 'absence' %in% presence_absence) %>%
  ungroup()

# Recode to 'yes' and 'no' for easier interpretation
biomass_df$absence_in_block[biomass_df$absence_in_block] = 'yes'
biomass_df$absence_in_block[biomass_df$absence_in_block != 'yes'] = 'no'

# 2.7 Set up outer cross-validation --------------------------------------------

set.seed(42)
outer_folds = group_vfold_cv(biomass_df, v = n_outer_folds, group = group_var)

print('Outer folds created:')
print(outer_folds)
cat('\n')

# ==============================================================================
# 3. PREPARE MODEL =============================================================
# ==============================================================================

# 3.1 Model specification ------------------------------------------------------

# Uses default parameters: 
# https://parsnip.tidymodels.org/reference/details_boost_tree_lightgbm.html
# Most important predictors to tune:
# https://lightgbm.readthedocs.io/en/latest/Parameters-Tuning.html
# https://github.com/microsoft/LightGBM/issues/695
# https://www.r-bloggers.com/2020/08/how-to-use-lightgbm-with-tidymodels/
lgbm = boost_tree( 
  tree_depth = tune(),
  min_n = tune()
) %>%
  set_mode("classification") %>%
  set_engine("lightgbm",
             num_leaves = tune())

# Uses default parameters: 
# https://parsnip.tidymodels.org/reference/details_rand_forest_ranger.html
rf = rand_forest( 
  mtry = tune()
) %>%
  set_mode("classification") %>%
  set_engine("ranger", 
             importance = "permutation",
             max.depth = tune())

# Uses default parameters: 
# https://parsnip.tidymodels.org/reference/details_svm_poly_kernlab.html
svmLinear = svm_linear( # Calls kernlab::ksvm() with kernel = "vanilladot"
  cost = tune()
) %>% 
  set_mode("classification") %>% # Default type is C-svc for classification
  set_engine("kernlab")

# Uses default parameters: 
# https://parsnip.tidymodels.org/reference/details_svm_poly_kernlab.html
svmPoly = svm_poly( # Calls kernlab::ksvm() with kernel = "polydot"
  cost = tune(),
  degree = tune(),
  scale_factor = tune()
) %>% 
  set_mode("classification") %>% # Default type is C-svc for classification
  set_engine("kernlab")

# Uses default parameters: 
# https://parsnip.tidymodels.org/reference/details_boost_tree_xgboost.html
xgb = boost_tree(
  tree_depth = tune(),
  min_n = tune(),
  loss_reduction = tune()
) %>% 
  set_mode("classification") %>%
  set_engine("xgboost")

# 3.2 Model formula ------------------------------------------------------------
# Must define formula manually rather than using 'update_role' because 'geometry' from sf is stored as list

form = as.formula(paste("presence_absence ~ ", paste(predictors, collapse =" + "), sep = ""))

print(paste0(c('The model formula is: ', form), collapse = ' '))
cat('\n')

# 3.3 Set up recipe ------------------------------------------------------------
# Recipe steps are applied to CV folds, the full data set is just used here to initialize recipe:
# https://community.rstudio.com/t/fit-resamples-prep-and-outcome-transform-behavior/83891

# Create recipe iterations
# Change feature selection threshold for each iteration
rec_list = list()
for(n in feat_thres_list){
  
  rec_n = recipe(formula = form, data = biomass_df) %>%
    step_nzv(all_numeric_predictors()) %>% # Remove predictors with near zero variance
    step_select_mrmr(all_predictors(), # MRMR feature selection
                     outcome = 'presence_absence',
                     threshold = n) %>% 
    step_novel(all_factor_predictors(), new_level = 'other') %>% # Allow recipe to identify factor levels not seen in the training data and treat them as a new level called "other"
    step_dummy(all_factor_predictors(), one_hot = TRUE) # One hot encode to keep all factor levels for variable importance purposes
  
  rec_list[[paste0('feat_thresh_', n * 100)]] = rec_n
  
}

print('Test feature selection recipe:')
print(rec_list[[1]])
cat('\n')

# 3.4 Final model grid ---------------------------------------------------------

grid_ctrl =
  control_grid(
    save_pred = TRUE,
    parallel_over = "everything",
    save_workflow = TRUE
  )

# ==============================================================================
# 4. LOOP THROUGH MODELS =======================================================
# ==============================================================================

for(mod_name in model_list){ # START MODEL TYPE LOOP
  
  # 4.1 Set up workflow set ----------------------------------------------------
  
  # Add data normalization to recipe, if necessary
  if(mod_name %in% models_normalize){
    
    # Preserve original recipe list
    rec_list_final = rec_list
    
    # Get recipe names
    rec_names = names(rec_list_final)
    
    # Add normalization step
    for(rec_name in rec_names){
      
      rec_current = rec_list_final[[rec_name]]
      rec_update = rec_current %>% step_normalize(all_predictors())
      rec_list_final[[rec_name]] = rec_update
      
    }
    
  }else{
    rec_list_final = rec_list
  }
  
  # Create workflow set
  wf = workflow_set(
    preproc = rec_list,
    models = list(get(mod_name))
  )

  print(paste0('The current model is: ', mod_name))
  cat('\n')
  
  print('The workflow set is: ')
  print(wf)
  cat('\n')
  
  # 4.2 Get outer fold split ---------------------------------------------------
  
  # Get training and testing data
  outer_fold = outer_folds[outer_fold_num,]
  outer_train = analysis(outer_fold$splits[[1]])
  outer_test = assessment(outer_fold$splits[[1]])
  
  print(paste0('Tuning outer fold ', outer_fold_num, ' of ', nrow(outer_folds)))
  cat('\n')
  
  # 4.3 Set up inner folds -----------------------------------------------------
  
  set.seed(11) # Set seed for reproducible results
  inner_folds = group_vfold_cv(outer_train, v = n_inner_folds, group = group_var, strata = strata_var)
  
  print(paste0(n_inner_folds, ' inner folds created:'))
  print(inner_folds)
  cat('\n')
  
  # 4.4 Tune over inner folds --------------------------------------------------

  # Adjust parameters if necessary
  if(mod_name == 'rf' | mod_name == 'rf_log'){
    
    # For random forest, manually update to initialize max.depth parameter
    param_info =
      rf %>%
      parameters() %>%
      update(max.depth = tree_depth())
    
  }else{
    
    # For all other models, use the default parameters
    param_info = NULL
    
  }

  # Tune over hyperparameters and number of features
  # Default function applied across the workflows is tune_grid():
  # https://workflowsets.tidymodels.org/reference/workflow_map.html
  # If no tuning grid is provided, a semi-random grid (via dials::grid_latin_hypercube()) is created:
  # https://tune.tidymodels.org/reference/tune_grid.html
  grid_results =
    wf %>%
    workflow_map(
      seed = 65,
      resamples = inner_folds,
      grid = grid_n,
      param_info = param_info,
      metrics = metrics,
      verbose = TRUE,
      control = grid_ctrl
    )
  
  # Get tidymodels model name
  tidy_mod_name = grid_results$info[[1]]$model
  
  # 4.5 Loop through results ---------------------------------------------------
  # Each row is a different model configuration with a different feature selection threshold
  # Get hyperparameters used in latin hypercube tuing
  # Get total number of features selected during feature selection
  
  # Initialize metadata dataframes
  parameter_grid_lookup = data.frame()
  n_feat_lookup = data.frame()
  
  for(j in 1:nrow(grid_results)){ # START FEATURE SELECTION THRESHOLD LOOP
    
    # Get hyperparameters used in latin hypercube tuning
    wf_row = grid_results[j,]
    wf_params = unique(wf_row$result[[1]]$.metrics[[1]] %>% dplyr::select(-c(.metric, .estimator, .estimate)))
    wf_params$wflow_id = wf_row$wflow_id
    
    # Add to metadata data frame
    parameter_grid_lookup = rbind(parameter_grid_lookup, wf_params)
    
    # Loop through inner cross-validation fold results
    # Get total number of features selected during feature selection
    for(k in 1:nrow(wf_row$result[[1]])){ # START INNER FOLD LOOP
      
      # Get data and feature selection threshold
      inner_fold_row = wf_row$result[[1]][k,]
      inner_fold_training_data = analysis(inner_fold_row$splits[[1]])
      inner_fold_feat_thres = gsub('feat_thresh_', '', wf_row$wflow_id)
      inner_fold_feat_thres = as.numeric(gsub(paste0('_', tidy_mod_name), '', inner_fold_feat_thres))/100
      
      # Prepare data according to recipe and count total number of features chosen
      inner_fold_n_feat = ncol(recipe(form, data = inner_fold_training_data) %>%
                                 step_nzv(all_numeric_predictors()) %>%
                                 step_select_mrmr(all_numeric_predictors(),
                                                  outcome = 'presence_absence',
                                                  threshold = inner_fold_feat_thres) %>%
                                 prep() %>%
                                 juice())
      
      # Tidy
      n_feat = data.frame(wflow_id = wf_row$wflow_id, id = inner_fold_row$id, n_feat = inner_fold_n_feat)
      
      # Add to metadata data frame
      n_feat_lookup = rbind(n_feat_lookup, n_feat)
      
    } # END INNER FOLD LOOP
    
  } # END FEATURE SELECTION THRESHOLD LOOP
  
  # 4.6 Aggregate accuracy results ---------------------------------------------
  
  # Aggregate accuracy results across all inner cv folds
  outer_fold_metrics = collect_metrics(grid_results, summarize = FALSE)
  
  # Add feature threshold column
  outer_fold_metrics$feat_thresh = gsub('feat_thresh_', '', outer_fold_metrics$wflow_id)
  outer_fold_metrics$feat_thresh = gsub(paste0('_', tidy_mod_name), '', outer_fold_metrics$feat_thresh)
  
  # Join number of features column
  outer_fold_metrics = left_join(outer_fold_metrics,
                                 n_feat_lookup,
                                 by = c('id', 'wflow_id'))
  
  # Add outer fold identifier
  outer_fold_metrics$outer_fold = outer_fold_num
  
  # Join parameter values
  outer_fold_metrics = left_join(outer_fold_metrics,
                                 parameter_grid_lookup,
                                 by = c('.config', 'wflow_id'))
  
  # 4.7 Get best model ---------------------------------------------------------
  
  # Get best model i.e. optimal # of features, best hyperparameters
  final_wf = fit_best(grid_results, metric = selection_metric)
  
  # Get best model parameters
  final_wf_metadata = rank_results(grid_results, rank_metric = selection_metric, select_best = TRUE)[1,]
  final_wf_id = final_wf_metadata$wflow_id
  final_feat_thresh = as.numeric(strsplit(final_wf_metadata$wflow_id, '_')[[1]][3])
  final_config = final_wf_metadata$.config
  final_params = parameter_grid_lookup[parameter_grid_lookup$wflow_id == final_wf_id & parameter_grid_lookup$.config == final_config,]
  final_params = final_params %>% dplyr::select(-c(.config, wflow_id))
  final_params$feat_thresh = final_feat_thresh
  
  # Get best model predictors
  final_predictors = data.frame(predictor = names(final_wf$pre$mold$predictors), outer_fold = outer_fold_num)
  
  # 4.8 Predict on outer fold training set -------------------------------------
  # Only used to test for overfitting
  
  # Apply best model to outer fold training set
  train_predictions = predict(final_wf, outer_train, type = "prob")
  
  # Aggregate training set predictions
  train_predictions_df =  st_drop_geometry(outer_train)
  train_predictions_df$prob_presence = data.frame(train_predictions)$.pred_presence
  train_predictions_df$prob_absence = data.frame(train_predictions)$.pred_absence
  train_predictions_df$outer_fold = outer_fold_num
  train_predictions_df$model = mod_name
  
  # Get class predictions from probabilities
  train_predictions_df = train_predictions_df %>%
    mutate(
      predicted = make_two_class_pred(
        estimate = prob_presence, 
        levels = levels(presence_absence), 
        threshold = .5
      ))
  
  # 4.9 Predict on outer fold test set -----------------------------------------
  # Used for final accuracy results
  
  # Apply best model to outer fold test set
  test_predictions = predict(final_wf, outer_test, type = "prob")
  
  # Aggregate test set predictions
  test_predictions_df =  st_drop_geometry(outer_test)
  test_predictions_df$prob_presence = data.frame(test_predictions)$.pred_presence
  test_predictions_df$prob_absence = data.frame(test_predictions)$.pred_absence
  test_predictions_df$outer_fold = outer_fold_num
  test_predictions_df$model = mod_name
  
  # Get class predictions from probabilities
  test_predictions_df = test_predictions_df %>%
    mutate(
      predicted = make_two_class_pred(
        estimate = prob_presence, 
        levels = levels(presence_absence), 
        threshold = .5
      ))
  
  # Attach optimal hyperparameters
  test_predictions_df = cbind(test_predictions_df, final_params)
  
  # 4.10 Aggregate training and test set accuracy results ----------------------
  # Used to test for overfitting
  
  compare_train = train_predictions_df
  compare_test = test_predictions_df
  
  # Aggregate training set metrics
  train_metrics = metrics(data = compare_train, truth = presence_absence, estimate = predicted, prob_presence)
  train_metrics$set = 'train'
  train_metrics$outer_fold = outer_fold_num
  train_metrics$model = mod_name
  
  # Aggregate test set metrics
  test_metrics = metrics(data = compare_test, truth = presence_absence, estimate = predicted, prob_presence)
  test_metrics$set = 'test'
  test_metrics$outer_fold = outer_fold_num
  test_metrics$model = mod_name
  
  # Combine
  compare_train_test = rbind(train_metrics, test_metrics)
  
  # Attach optimal hyperparameters
  compare_train_test = cbind(compare_train_test, final_params)
  
  # Specify classification threshold
  compare_train_test$class_threshold = 0.5
  
  # 4.11 Calculate accuracy for different classification thresholds ------------
  
  # Set up F measure metric with beta parameter
  f_meas_beta2 = metric_tweak("f_meas_beta2", f_meas, beta = 2)
  
  # Calculate threshold performances
  thresholds = test_predictions_df %>%
    threshold_perf(truth = presence_absence, 
                   estimate = prob_presence, 
                   thresholds = seq(0, 1, by = 0.05),
                   metrics = metric_set(yardstick::accuracy, yardstick::bal_accuracy, yardstick::kap, yardstick::f_meas, f_meas_beta2, yardstick::j_index, yardstick::mcc, yardstick::sens, yardstick::spec, yardstick::precision))
  
  # Tidy threshold data
  thresholds = thresholds %>% rename(threshold = .threshold, metric = .metric, estimator = .estimator, estimate = .estimate)
  thresholds$outer_fold = outer_fold_num
  thresholds$model = mod_name
  
  # Attach optimal hyperparameters
  thresholds = cbind(thresholds, final_params)
  
  # 4.12 Add accuracy and final predictions to model data frame ----------------
  
  # Convert predictions to factor to facilitate rbind
  test_predictions_df$predicted = as.factor(test_predictions_df$predicted)
  
  # Remove predicted class from final test predictions, this will be determined later using optimal threshold
  test_predictions_df = test_predictions_df %>% dplyr::select(-predicted)

  # 4.13 Save ------------------------------------------------------------------
  
  # Save metrics for current model
  write.csv(outer_fold_metrics, paste0(out_dir, '/accuracy_metrics/accuracy_metrics_', ds_type, '_', mod_name, '_', outer_fold_num, '.csv'), row.names = FALSE)
  
  # Save final predictors for current model
  write.csv(final_predictors, paste0(out_dir, '/predictors/final_predictors_', ds_type, '_', mod_name, '_', outer_fold_num, '.csv'), row.names = FALSE)
  
  # Save final predictions for current model (on outer fold test set)
  write.csv(test_predictions_df, paste0(out_dir, '/predictions/final_predictions_', ds_type, '_', mod_name, '_', outer_fold_num, '.csv'), row.names = FALSE)
  
  # Save final train/test comparison for current model
  write.csv(compare_train_test, paste0(out_dir, '/train_test/compare_train_test_', ds_type, '_', mod_name, '_', outer_fold_num, '.csv'), row.names = FALSE)
  
  # Save final threshold comparison for current model
  write.csv(thresholds, paste0(out_dir, '/thresholds/classification_thresholds_', ds_type, '_', mod_name, '_', outer_fold_num, '.csv'), row.names = FALSE)

} # END MODEL TYPE LOOP
