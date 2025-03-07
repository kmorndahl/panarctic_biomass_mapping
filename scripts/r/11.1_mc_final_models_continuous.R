################################################################################
################################################################################

# DESCRIPTION:

# This script runs final regression model Monte Carlo iterations (total and woody):
# - data grouped by spatial blocks
# - data stratified by presence/absence to prevent class imbalance
# - cross-validation with 10 folds
# - random forest model, tuned mtry, max.depth
# - MRMR feature selection
# - Unvegetated data used for training

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

# REFERENCES:

# - Cross-validation best practices: https://hastie.su.domains/MOOC-Slides/cv_boot.pdf
# - Nested CV example: https://www.tidymodels.org/learn/work/nested-resampling/
# - Recommended pre-processing: https://www.tmwr.org/pre-proc-table.html
# - All recipe step functions: https://recipes.tidymodels.org/reference/index.html
# - Hurdle model example: https://geoffruddock.com/building-a-hurdle-regression-estimator-in-scikit-learn/
# - Workflow sets: https://www.tmwr.org/workflow-sets.html

# NOTES:
# - ranger version 0.15.4 or higher reqiured to use 'node.stats' parameter
#   - https://imbs-hl.github.io/ranger/news/index.html

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

# library(tidyverse)
library(dplyr)
library(tidyr)
library(purrr)
library(classInt)

# library(tidymodels)
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

library(vip)

source('scripts/r/000.0_yardstick_bias.R')
source('scripts/r/000.1_yardstick_rmspe.R')

# 1.2 Parameters ---------------------------------------------------------------

# Bash
params = commandArgs(trailingOnly = TRUE)
ds_type = params[1] # Get dataset type from bash script
i = c(params[2]) # Get monte-carlo iteration from bash script
 
# R
in_version = 'v20240508'
out_version = 'v20240514'
mod_name = 'rf'
presence_thresh = 10 # Presence threshold
metrics = metric_set(yardstick::rmse, yardstick::mae, yardstick::mape, yardstick::rsq, rmspe, bias) # Evaluation metrics
grid_n = 20 # Hyperparameter tuning grid size
n_folds = 10
feat_thres_list = c(0.50, 0.75, 0.90, 0.95) # Feature threshold options
selection_metric = 'rmse' # Metric for choosing best parameters - choose 'rmse' or 'mae'
group_var = 'spatial_block_id' # CV grouping variable
strata_var = 'biomass_bin' # CV stratification variable
plot_area_threshold = 0.05 # Site level plot area threshold
ccdc_rmse_threshold = 2500 # Site level CCDC RMSE threshold

# Directories
in_dir = 'output/06_model_ready_data/plot_level_mc/'
out_dir = 'output/11_mc_final_models/continuous/'

# 1.3 Functions ----------------------------------------------------------------

find_mode = function(x) {
  u = unique(na.omit(x))
  tab = tabulate(match(x, u))
  u[tab == max(tab)][1] # If there are two modes, grab the first one
}

std.error = function(x) sd(na.omit(x))/sqrt(length(na.omit(x)))

# 1.4 Read in data -------------------------------------------------------------

# Read in PLOT-LEVEL synthesis dataset
biomass_df = read.csv(paste0(in_dir, 'ds_', ds_type, '_plots_filtered_', in_version, '.csv'))

# 1.5 Report -------------------------------------------------------------------

cat('\n')
print(paste0('In version: ', in_dir))
print(paste0('Out version: ', out_dir))
print(paste0('Data set type: ', ds_type))
print(paste0('Model: ', mod_name))
print(paste0('Presence threshold: ', presence_thresh))
print(paste0('Hyperparameter grid size: ', grid_n))
print(paste0('Number of folds: ', n_folds))
print(paste0('Feature threshold options: ', paste(feat_thres_list, collapse = ', ')))
print(paste0('Selection metric: ', selection_metric))
print(paste0('Grouping variable: ', group_var))
print(paste0('Stratification variable: ', strata_var))
cat('\n')

################################################################################
################################################################################

# ==============================================================================
# 2. TIDY DATA =================================================================
# ==============================================================================

# 2.1 Remove unnecessary columns -----------------------------------------------

# Specify columns to remove
remove_cols = grep("biomass_dry_weight_g|rep_|pixel_num|distance|n_plots|n_pixels", names(biomass_df), value = TRUE)

# Specify additional columns to remove
# Remove spectral summaries - we need to recalculate these using permuted data
remove_cols_spectral_summaries = names(biomass_df)[grepl("_change|_annualMean|_annualMedian|_annualRange", names(biomass_df))]
remove_cols_spectral_summaries = remove_cols_spectral_summaries[remove_cols_spectral_summaries != "ccdc_changeProb"] # Retain ccdc_changeProb

# Remove unneeded columns
biomass_df = biomass_df %>% select(-any_of(remove_cols))
biomass_df = biomass_df %>% select(-any_of(remove_cols_spectral_summaries))

# 2.2 Specify column handling --------------------------------------------------

# Specify columns to average
average_cols = grep("ccdc_|coefs_|doy_|spectral_|texture_|climate_|topo_|soil_|trees_cover|permafrost|biomass_density_gm2|plot_area_m2_mean|latitude|longitude", names(biomass_df), value = TRUE)
average_cols = average_cols[average_cols != "ccdc_segmentType"] # Remove segment type, it should not be averaged

# Convert columns to numeric if they are to be averaged
# There will be some NAs introduced for: ccdc_NextSegmentStart, ccdc_NextSegmentEnd, ccdc_previousSegmentStart, ccdc_previousSegmentEnd
# That is okay, this is where these are 'no segment'
biomass_df = biomass_df %>% mutate_at(average_cols, as.numeric)

# 2.3 Drop NAs -----------------------------------------------------------------

# Specify columns where NAs can be dropped
drop_na_cols = grep("soil_|world_terrestrial_ecosystems", names(biomass_df), value = TRUE)

# A few NAs for soil and world_terrestrial_ecosystems
# These are nonveg points in ice or water, we can safely remove them
biomass_df = biomass_df %>% drop_na(any_of(drop_na_cols))

################################################################################
################################################################################

# ==============================================================================
# 3. SPECIFY PREDICTORS ========================================================
# ==============================================================================

# Initialize full list of predictors
predictors = grep("permafrost|spectral_|texture_|topo_|trees_|zone_3|ecoregion|world_terrestrial_ecosystems|latitude|longitude", names(biomass_df), value = TRUE)
predictors = predictors[!grepl("_se|_sd|harvest_", predictors)] # Remove standard error, standard deviation and harvest_doy

# Get CCDC RMSE
ccdc_metadata = st_drop_geometry(biomass_df) %>%
  select(matches("_rmse") & !matches("_se|_avg"))

################################################################################
################################################################################

# ==============================================================================
# 4. SAMPLE WITHIN SITE ========================================================
# ==============================================================================

# ******************************************************************************
# Sample plots/pixels within a site with replacement
# Accounts for sampling variability of plots across site
# Later, average or take mode of plots/pixels within site to get site-level data
# ******************************************************************************

set.seed(i)

# Aggregate data
biomass_df = biomass_df %>%
  group_by(site_code) %>% # Group by site
  slice_sample(prop = 1, replace = TRUE) # Sample with replacement within site

################################################################################
################################################################################

# ==============================================================================
# 5. PERMUTE LANDSAT PREDICTOR DATA ============================================
# ==============================================================================

# 5.1 Get spectral predictors --------------------------------------------------

bands_si = predictors[grepl('spectral_', predictors)]
bands_si = gsub('spectral_', '', bands_si)
bands_si = unique(gsub('_startSnowfree|_earlySummer|_peakSummer|_lateSummer|_endSnowfree', '', bands_si))

# 5.2 Loop bands/indices -------------------------------------------------------

for(b in bands_si){ # START BAND/INDEX LOOP
  
  # Get predictors for current band
  predictors_band_si = predictors[grepl(b, predictors)]
  
  # 5.3 CCDC model uncertainty -------------------------------------------------
  
  # ****************************************************************************
  # Permute spectral bands and indices using CCDC RMSE
  # Accounts for uncertainty in CCDC models
  # ****************************************************************************
  
  set.seed(i)
  
  # Generate n plot x n band random draws from standard normal distribution - each band/season/observation gets a random number
  ccdc_permutation = data.frame(matrix(rnorm(n = length(predictors_band_si) * nrow(biomass_df), mean = 0, sd = 1), ncol = length(predictors_band_si)))
  names(ccdc_permutation) = predictors_band_si
  
  # Get band CCDC RMSE
  ccdc_rmse = ccdc_metadata %>% select(matches(b))
  
  # Permute
  for(predictor in predictors_band_si){ # START PREDICTOR LOOP
    
    # Permute values
    permuted_values = biomass_df[predictor] + (ccdc_rmse * ccdc_permutation[predictor])
    
    # Clamp values so they are all still valid
    if(b %in%  c("NIR", "SWIR1", "SWIR2", "blue", "green", "red")){
      
      # If it is a spectral band, clamp to (0 , 10,000)
      permuted_values[permuted_values < 0] = 0
      permuted_values[permuted_values > 10000] = 10000
      
    }else if(b %in%  c("EVI2b", "NBR", "NDMI", "NDVI", "NDWI")){
      
      # If it is a spectral band, clamp to (-10,000 , 10,000)
      permuted_values[permuted_values < -10000] = -10000
      permuted_values[permuted_values > 10000] = 10000
      
    }else{
      print('Band not identified as band or index, no clamping applied')
    }
    
    # Replace original values with permuted values
    biomass_df[predictor] = permuted_values
    
  } # END PREDICTOR LOOP
  
  # 5.4 Calculate spectral summaries -------------------------------------------
  
  # Get spectral band names
  start_snowfree = paste0('spectral_', b, '_startSnowfree')
  early_summer = paste0('spectral_', b, '_earlySummer')
  peak_summer = paste0('spectral_', b, '_peakSummer')
  late_summer = paste0('spectral_', b, '_lateSummer')
  end_snowfree = paste0('spectral_', b, '_endSnowfree')
  
  # Get seasonal DOY names
  start_snowfree_doy = 'doy_startSnowfree'
  early_summer_doy = 'doy_earlySummer'
  peak_summer_doy = 'doy_peakSummer'
  late_summer_doy = 'doy_lateSummer'
  end_snowfree_doy = 'doy_endSnowfree'
  
  # Assign output band names
  annualMean = paste0('spectral_', b, '_annualMean')
  annualMedian = paste0('spectral_', b, '_annualMedian')
  annualMin = paste0('spectral_', b, '_annualMin')
  annualMax = paste0('spectral_', b, '_annualMax')
  annualRange = paste0('spectral_', b, '_annualRange')
  changeSSES = paste0('spectral_', b, '_changeSSES')
  changeESPS = paste0('spectral_', b, '_changeESPS')
  changePSLS = paste0('spectral_', b, '_changePSLS')
  changeLSES = paste0('spectral_', b, '_changeLSES')
  
  # Calculate annual summaries and seasonal rates of change
  # NOTE: take absolute value of day change for rates
  #   - for some barrens, DOYs are inconsistent i.e. sometimes peak summer is actually a day or so later than late summer
  #   - shouldn't be a long term problem, usually only a day difference, likely due to prevalence of snow/ice, lack of data
  biomass_df = biomass_df %>%
    rowwise() %>%
    mutate(
      !!annualMedian := median(c_across(all_of(c(start_snowfree, early_summer, peak_summer, late_summer, end_snowfree)))),
      !!annualMean := mean(c_across(all_of(c(start_snowfree, early_summer, peak_summer, late_summer, end_snowfree)))),
      !!annualMin := min(c_across(all_of(c(start_snowfree, early_summer, peak_summer, late_summer, end_snowfree)))),
      !!annualMax := max(c_across(all_of(c(start_snowfree, early_summer, peak_summer, late_summer, end_snowfree)))),
      !!annualRange := get(annualMax) - get(annualMin),
      !!changeSSES := (get(early_summer)-get(start_snowfree))/abs(get(early_summer_doy)-get(start_snowfree_doy)),
      !!changeESPS := (get(peak_summer)-get(early_summer))/abs(get(peak_summer_doy)-get(early_summer_doy)),
      !!changePSLS := (get(late_summer)-get(peak_summer))/abs(get(late_summer_doy)-get(peak_summer_doy)),
      !!changeLSES := (get(end_snowfree)-get(late_summer))/abs(get(end_snowfree_doy)-get(late_summer_doy))
    ) %>%
    select(-c(annualMin, annualMax))
  
  # Replace NaN and Inf with zeros, this is where slopes have 0 denominators
  biomass_df = biomass_df %>% mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))
  biomass_df = biomass_df %>% mutate(across(everything(), ~replace(.x, is.infinite(.x), 0)))
  
} # END BAND/INDEX LOOP

################################################################################
################################################################################

# ==============================================================================
# 6. AGGREGATE TO SITE LEVEL ===================================================
# ==============================================================================

# Specify factors
biomass_df = biomass_df %>% mutate_at(c('zone_3', 'ecoregion', 'land_cover', 'world_terrestrial_ecosystems'), as.factor)

# Aggregate
biomass_df = biomass_df %>%
  group_by(site_code) %>% # Group by site
  summarise(
    across(all_of(average_cols), mean, na.rm = TRUE), # For numeric columns, get mean
    across(-average_cols, find_mode) # For non-numeric columns, get mode
  )

################################################################################
################################################################################

# ==============================================================================
# 7. FILTER SITE LEVEL DATA ====================================================
# ==============================================================================

# Calculate average CCDC RMSE across all bands/indices
biomass_df = biomass_df %>% mutate(biomass_df, ccdc_rmse_avg = rowMeans(select(biomass_df, ends_with("_rmse")), na.rm = TRUE))

# Plot area
biomass_df = biomass_df[biomass_df$plot_area_m2_mean >= plot_area_threshold,]

# CCDC RMSE
biomass_df = biomass_df[biomass_df$ccdc_rmse_avg < ccdc_rmse_threshold,]

################################################################################
################################################################################

# ==============================================================================
# 8. PREPARE DATA FOR MODELING =================================================
# ==============================================================================

# Get lat/long
latlong = biomass_df %>% dplyr::select(any_of(c("longitude", "latitude")))

# Convert data to sf
biomass_df = sf::st_as_sf(biomass_df, coords = c("longitude", "latitude"), crs = 4326)
biomass_df = cbind(biomass_df, latlong)

# Transform biomass data
biomass_df = st_transform(biomass_df, 3571)

################################################################################
################################################################################

# ==============================================================================
# 9. SET UP SPATIAL CROSS-VALIDATION ===========================================
# ==============================================================================

# 9.1 Create spatial blocks ----------------------------------------------------

set.seed(1908) # Set seed for reproducible results
block_folds = spatialsample::spatial_block_cv(biomass_df, v = NULL, cellsize = 100*1000) # Cellsize units should be in m: 100 km times 1000 m in a km

print('Spatial block folds created:')
print(block_folds)
cat('\n')

# 9.2 Assign spatial block IDs to data -----------------------------------------

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

################################################################################
################################################################################

# ==============================================================================
# 10. SAMPLE OBSERVATIONS WITHIN DATASET =======================================
# ==============================================================================

# 10.1 Spatial variability across sites ----------------------------------------

# ******************************************************************************
# Sample sites within dataset with replacement
# Accounts for sampling variability of sites across study area
# ******************************************************************************

set.seed(i)

# Sample with replacement across sites
# Accounts for: sampling error across sites
biomass_df = dplyr::sample_n(biomass_df, nrow(biomass_df), replace = TRUE)

# 10.2 Grab new predictor list -------------------------------------------------

# Initialize full list of predictors
predictors = grep("permafrost|spectral_|texture_|topo_|trees_|zone_3|ecoregion|world_terrestrial_ecosystems|latitude|longitude", names(biomass_df), value = TRUE)
predictors = predictors[!grepl("_se|_sd|harvest_", predictors)] # Remove standard error, standard deviation and harvest_doy

################################################################################
################################################################################

# ==============================================================================
# 11. BIN RESPONSE VARIABLE =========================================
# ==============================================================================
# NOTE: observations are assigned 'high' biomass if ANY observation in their spatial block is >= the threshold
# - This is necessary so that the CV folds can be both grouped by spatial block and stratified by biomass

# Calculate intervals using presence only data
biomass_breaks = classIntervals(biomass_df[biomass_df$biomass_density_gm2>0,]$biomass_density_gm2, n=3, style="fisher") # 'low', 'medium' and 'high' -- label 'medium' and 'high' as 'high' to ensure enough blocks to distribute among cross-validation folds

# Get the cut-off value and round down to nearest 100th
biomass_bin_cutoff = floor(biomass_breaks$brks[2]/100)*100

# Biomass >= biomass_bin_cutoff g/m2 is 'high'
biomass_df$biomass_bin = 'low'
biomass_df$biomass_bin[biomass_df$biomass_density_gm2 >= biomass_bin_cutoff] = 'high'

# Assign 'high' biomass to a spatial block if any observation within that block has high biomass
biomass_df = biomass_df %>%
  group_by(spatial_block_id) %>%
  mutate(biomass_bin = 'high' %in% biomass_bin) %>%
  ungroup()

# Recode to 'low' and 'high' for easier interpretation
biomass_df$biomass_bin[biomass_df$biomass_bin] = 'high'
biomass_df$biomass_bin[biomass_df$biomass_bin != 'high'] = 'low'

################################################################################
################################################################################

# ==============================================================================
# 12. MODEL SET UP =============================================================
# ==============================================================================

rf = rand_forest( # Tuning parameters: https://parsnip.tidymodels.org/reference/details_rand_forest_ranger.html
  mtry = tune()
) %>%
  set_mode("regression") %>%
  set_engine("ranger", 
             importance = "permutation", 
             node.stats = TRUE,
             max.depth = tune())

# Formula
# Must define formula manually rather than using 'update_role' because 'geometry' from sf is stored as list
form = as.formula(paste("biomass_density_gm2 ~ ", paste(predictors, collapse =" + "), sep = ""))
print(paste0(c('The model formula is: ', form), collapse = ' '))
cat('\n')

################################################################################
################################################################################

# ==============================================================================
# 13. RECIPE SET UP ============================================================
# ==============================================================================

# Note: recipe steps are applied to CV folds, the full data set is just used here to initialize recipe
#   - https://community.rstudio.com/t/fit-resamples-prep-and-outcome-transform-behavior/83891

# Create recipe iterations changing feature selection number
rec_list = list()
for(n in feat_thres_list){
  
  rec_n = recipe(formula = form, data = biomass_df) %>%
    step_nzv(all_numeric_predictors()) %>%
    step_select_mrmr(all_predictors(),
                     outcome = 'biomass_density_gm2',
                     threshold = n) %>% 
    step_novel(all_factor_predictors(), new_level = 'other') %>% # Allow recipe to identify factor levels not seen in the training data and treat them as a new level called "other"
    step_dummy(all_factor_predictors(), one_hot = TRUE) # One hot encode to keep all factor levels for variable importance purposes

  rec_list[[paste0('feat_thresh_', n * 100)]] = rec_n
  
}

print('Test feature selection recipe:')
print(rec_list[[1]])
cat('\n')

################################################################################
################################################################################

# ==============================================================================
# 14. WORKFLOW SET UP ==========================================================
# ==============================================================================

grid_ctrl =
  control_grid(
    save_pred = TRUE,
    parallel_over = "everything",
    save_workflow = TRUE
  )

# Set up workflow set

wf = workflow_set(
  preproc = rec_list,
  models = list(get(mod_name))
)
response_var = 'biomass_density_gm2'

print('The workflow set is: ')
print(wf)
cat('\n')

################################################################################
################################################################################

# ==============================================================================
# 15. CROSS VALIDATION SET UP ==================================================
# ==============================================================================

set.seed(42) # Set seed for reproducible results
folds = group_vfold_cv(biomass_df, v = n_folds, group = group_var, strata = strata_var)

print(paste0(n_folds, ' folds created:'))
print(folds)
cat('\n')

################################################################################
################################################################################

# ==============================================================================
# 16. TUNE OVER FOLDS ==========================================================
# ==============================================================================

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
# Default function applied across the workflows is tune_grid()
# https://workflowsets.tidymodels.org/reference/workflow_map.html
# If no tuning grid is provided, a semi-random grid (via dials::grid_latin_hypercube()) is created
# https://tune.tidymodels.org/reference/tune_grid.html
grid_results =
  wf %>%
    workflow_map(
      seed = 65,
      resamples = folds,
      grid = grid_n,
      param_info = param_info,
      metrics = metrics,
      verbose = TRUE,
      control = grid_ctrl
    )

# Get tidymodels model name
tidy_mod_name = grid_results$info[[1]]$model

################################################################################
################################################################################

# ==============================================================================
# 17. GET PARAMETERS AND ASSOCIATED WORKFLOW IDS AND CONFIGURATIONS ============
# ==============================================================================

parameter_grid_lookup = data.frame()
n_feat_lookup = data.frame()

for(j in 1:nrow(grid_results)){
  
  wf_row = grid_results[j,]
  wf_params = unique(wf_row$result[[1]]$.metrics[[1]] %>% dplyr::select(-c(.metric, .estimator, .estimate)))
  wf_params$wflow_id = wf_row$wflow_id
  
  parameter_grid_lookup = rbind(parameter_grid_lookup, wf_params)
  
  if(mod_name != 'linearBaseline_log'){
    
    for(k in 1:nrow(wf_row$result[[1]])){
      
      inner_fold_row = wf_row$result[[1]][k,]
      inner_fold_training_data = analysis(inner_fold_row$splits[[1]])
      inner_fold_feat_thres = gsub('feat_thresh_', '', wf_row$wflow_id)
      inner_fold_feat_thres = as.numeric(gsub(paste0('_', tidy_mod_name), '', inner_fold_feat_thres))/100
      inner_fold_n_feat = ncol(recipe(form, data = inner_fold_training_data) %>%
                                 step_nzv(all_numeric_predictors()) %>%
                                 step_select_mrmr(all_numeric_predictors(),
                                                  outcome = 'biomass_density_gm2',
                                                  threshold = inner_fold_feat_thres) %>%
                                 prep() %>%
                                 juice())
      
      n_feat = data.frame(wflow_id = wf_row$wflow_id, id = inner_fold_row$id, n_feat = inner_fold_n_feat)
      
      n_feat_lookup = rbind(n_feat_lookup, n_feat)
      
    }
    
  }
  
}

################################################################################
################################################################################

# ==============================================================================
# 18. AGGREGATE METRICS ========================================================
# ==============================================================================

# Aggregate metrics across all inner cv folds
fold_metrics = collect_metrics(grid_results, summarize = FALSE)

# Add feature threshold column
fold_metrics$feat_thresh = gsub('feat_thresh_', '', fold_metrics$wflow_id)
fold_metrics$feat_thresh = gsub(paste0('_', tidy_mod_name), '', fold_metrics$feat_thresh)

# Join number of features column
fold_metrics = left_join(fold_metrics,
                         n_feat_lookup,
                         by = c('id', 'wflow_id'))

# Tidy feature selection columns
if(grepl('_log', mod_name)){ # Test if model has log transformed response
  fold_metrics$log = 'yes'
  fold_metrics$n_feat = gsub('log_', '', fold_metrics$n_feat)
}else{
  fold_metrics$log = 'no'
}

# Join parameter values
fold_metrics = left_join(fold_metrics,
                         parameter_grid_lookup,
                         by = c('.config', 'wflow_id'))

################################################################################
################################################################################

# ==============================================================================
# 19. GET BEST MODEL ===========================================================
# ==============================================================================

# Get best model
# Fit using best workflow i.e. # of features, and best hyperparameters
set.seed(0)
final_wf = fit_best(grid_results, metric = selection_metric, verbose = 1)

# Get best model parameters
final_wf_metadata = rank_results(grid_results, rank_metric = selection_metric, select_best = TRUE)[1,]
final_wf_id = final_wf_metadata$wflow_id
final_feat = as.numeric(strsplit(final_wf_metadata$wflow_id, '_')[[1]][3])
final_config = final_wf_metadata$.config
final_params = parameter_grid_lookup[parameter_grid_lookup$wflow_id == final_wf_id & parameter_grid_lookup$.config == final_config,]
final_params = final_params %>% dplyr::select(-c(.config, wflow_id))
final_params$feat = final_feat
final_params$model = mod_name

################################################################################
################################################################################

# ==============================================================================
# 20. PREDICT ON FULL DATA SET =================================================
# ==============================================================================

# Apply best model to full data set
predictions = predict(final_wf, biomass_df)

# Aggregate training set predictions
predictions_df =  st_drop_geometry(biomass_df)
predictions_df$predicted = data.frame(predictions)$.pred

################################################################################
################################################################################

# ==============================================================================
# 21. GET PREDICTORS AND MODEL =================================================
# ==============================================================================

# Get predictors
final_predictors = data.frame(predictor = names(final_wf$pre$mold$predictors))

# Get parsnip model
parsnip_mod = final_wf %>% extract_fit_parsnip()

# Get model in original class (i.e. from engine used to fit model)
engine_mod = parsnip_mod$fit

################################################################################
################################################################################

# ==============================================================================
# 22. JUICE FINAL DATA =========================================================
# ==============================================================================

# Final recipe for prepping data
final_rec = recipe(formula = update(form,    ~ . + site_code), data = biomass_df) %>%
  update_role(site_code, new_role = "id variable") %>% # Save site code for identification
  step_nzv(all_numeric_predictors()) %>%
  step_select_mrmr(all_predictors(),
                   outcome = 'biomass_density_gm2',
                   threshold = final_params$feat / 100) %>% 
  step_novel(all_factor_predictors(), new_level = 'other') %>% # Allow recipe to identify factor levels not seen in the training data and treat them as a new level called "other"
  step_dummy(all_factor_predictors(), one_hot = TRUE) # One hot encode to keep all factor levels for variable importance purposes

# Final data
final_data_juiced = data.frame(final_rec %>% prep() %>% juice())

################################################################################
################################################################################

# ==============================================================================
# 23. VARIABLE IMPORTANCE ======================================================
# ==============================================================================
# https://koalaverse.github.io/vip/articles/vip.html#model-agnostic-vi
# https://github.com/koalaverse/vip/issues/131

# Function to get probability predictions from ranger random forest model
pred_fun = function(object, newdata) {
  predict(object, newdata, type = "response")$predictions
}

# Model specific importance, for RF we specified 'permutation'
vi_perm = parsnip_mod %>% 
  vi() 

# SHAP variable importance
vi_shap = engine_mod %>% 
  vi(method = 'shap',
     train = final_data_juiced,
     pred_wrapper = pred_fun)

################################################################################
################################################################################

# ==============================================================================
# 24. SAVE ALL =================================================================
# ==============================================================================

write.csv(final_predictors, paste0(out_dir, '/predictors/final_predictors_', ds_type, '_continuous_', out_version, '_', i,  '.csv'), row.names = FALSE)
write.csv(predictions_df, paste0(out_dir, '/predictions/final_predictions_', ds_type, '_continuous_', out_version, '_', i, '.csv'), row.names = FALSE)
write.csv(final_params, paste0(out_dir, '/parameters/final_parameters_', ds_type, '_continuous_', out_version, '_', i,  '.csv'), row.names = FALSE)
saveRDS(engine_mod, paste0(out_dir, '/models/', ds_type, '_continuous_', out_version, '_', i,  '.rds'))
write.csv(final_data_juiced, paste0(out_dir, '/final_data/final_data_', ds_type, '_continuous_', out_version, '_', i,  '.csv'), row.names = FALSE)
write.csv(vi_perm, paste0(out_dir, '/vip/vip_permutation_', ds_type,  '_continuous_', out_version, '_', i,  '.csv'), row.names = FALSE)
write.csv(vi_shap, paste0(out_dir, '/vip/vip_shap_', ds_type,  '_continuous_', out_version, '_', i,  '.csv'), row.names = FALSE)
