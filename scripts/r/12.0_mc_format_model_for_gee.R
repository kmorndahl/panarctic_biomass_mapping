################################################################################
################################################################################

# DESCRIPTION:

# This script converts random forest models into format readable by Google Earth Engine

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
library(ranger)
library(randomForest)
library(caret)
library(tidymodels)
library(gtools)
source('scripts/000.2_reprtree_update.R')

# 1.2 Parameters ---------------------------------------------------------------

options(scipen=999) # Ensure numbers are not represented in scientific notation
version = 'v20240514'
ds_types = c('total', 'woody') # Choose from 'total', 'woody'
response_types = c('continuous', 'binary') # Choose from 'continuous', 'binary'
model_type = 'ranger' # Specify R package used to create model - 'randomForest' or 'ranger'
model_fit_package = 'base' # Choose 'base' for models fit directly with 'randomForest' or 'ranger', choose 'caret' for models fit using 'caret', choose 'tidymodels' for models fit using tidymodels

# ==============================================================================
# 2. CONVERT MODELS ============================================================
# ==============================================================================

for(ds_type in ds_types){ # START DATASET TYPE LOOP
  
  print(paste0('===== The dataset type is: ', ds_type, ' ====='))
  cat('\n')
  
  for(response_type in response_types){ # START RESPONSE TYPE LOOP
    
    print(paste0('===== The response type is: ', response_type, ' ====='))
    cat('\n')
    
    # 2.1 Set up directories and parameters ------------------------------------
    
    in_dir = paste0('output/11_mc_final_models/', response_type, '/models/')
    out_dir = paste0('output/12_mc_trees/', response_type, '/')
    response_type_conversion = if(response_type == 'continuous') 'regression' else 'probability' 
    
    # 2.2 Get models -------------------------------------------------------------
    
    mod_list = list.files(in_dir, paste0(ds_type, '_', response_type, '_', version, '.*.rds'), full.names = TRUE)
    mod_list = mixedsort(mod_list) # Sort mc numbers properly
    mod_list = lapply(mod_list, readRDS)
    
    for(mc in 1:100){ # START MONTE CARLO ITERATION LOOP
      
      print(paste0('===== The Monte Carlo iteration is: ', mc, ' ====='))
      cat('\n')
      
      # 2.3 Get Monte Carlo model ------------------------------------------------
      
      mod = mod_list[[mc]]
  
      # 2.4 Convert forest -------------------------------------------------------
      
      out_name = paste0(out_dir, ds_type, '_', response_type, '_gee_', version, '_', mc, '.txt')
      mod = prep.mod(mod, model_type, response_type_conversion, model_fit_package)
      convert.forest(mod, out_name)
  
    } # END MONTE CARLO ITERATION LOOP

  } # END RESPONSE TYPE LOOP
    
} # END DATASET TYPE LOOP