################################################################################
################################################################################

# DESCRIPTION:

# This script converts compares predictions made in R and those made in 
# Google Earth Engine to ensure they are identical

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

# NOTES:

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

require(tidyverse)

# 1.2 Parameters ---------------------------------------------------------------

ds_type = 'woody'
response_type = 'continuous'
mod_type = 'rf'
version = 'v20240514'
mc = 1

# Directories
r_dir = paste0('output/11_mc_final_models/', response_type, '/predictions/')
gee_dir = paste0('output/13_mc_gee_compare/')

# 1.3 Read in data -------------------------------------------------------------

predictions_gee = read.csv(paste0(gee_dir, 'GEE_predictions_', response_type, '_', mod_type, '_', ds_type, '_', version, '_', mc, '.csv'))
predictions_r = read.csv(paste0(r_dir, 'final_predictions_', ds_type, '_', response_type, '_', version, '_', mc, '.csv'))

################################################################################
################################################################################

# ==============================================================================
# 2. TIDY DATA =================================================================
# ==============================================================================

# 2.1 Remove unnecessary columns -----------------------------------------------

predictions_gee = predictions_gee %>% select(-c('system.index', '.geo'))
predictions_r = predictions_r %>% select(all_of(names(predictions_gee)))
predictions_gee = predictions_gee %>% rename_at(vars(-contains( "site_code" )), ~paste0(., "_gee"))

# 2.2 Remove duplicates --------------------------------------------------------
# These result from sampling with replacement
# Predictions are identical across duplicated sites

predictions_gee = predictions_gee %>% distinct(site_code, .keep_all = TRUE)
predictions_r = predictions_r %>% distinct(site_code, .keep_all = TRUE)

# 2.3 Join datasets ------------------------------------------------------------
predictions = left_join(predictions_r, predictions_gee, by = 'site_code')

################################################################################
################################################################################

# ==============================================================================
# 3. PLOT ======================================================================
# ==============================================================================

if(response_type == 'continuous'){
  
  plt = ggplot(predictions, aes(y = predicted, x = predicted_gee))+
    geom_point(size = 2)+
    geom_smooth(method = 'lm', col = 'black', size = 0.5)+
    geom_abline(slope = 1, intercept = 0, lty = 2, size = 0.5)+
    labs(y = expression(paste('Total Biomass (g  ', m^-2, ') Predicted (R)')), x = expression(paste('Total Biomass (g  ', m^-2, ') Predicted (GEE)')))+
    theme_minimal(base_size = 30)
  plt

}else{
  
  plt = ggplot(predictions, aes(y = prob_presence, x = prob_presence_gee))+
    geom_point(size = 2)+
    geom_smooth(method = 'lm', col = 'black', size = 0.5)+
    geom_abline(slope = 1, intercept = 0, lty = 2, size = 0.5)+
    labs(y = 'Probability of Presence (R)', x = 'Probability of Presence (GEE)')+
    theme_minimal(base_size = 30)
  plt
  
}

################################################################################
################################################################################

# ==============================================================================
# 4. SAVE ======================================================================
# ==============================================================================

ggsave(
  paste0(gee_dir, 'r_gee_compare_', ds_type, '_', response_type, '_mc', mc, '.png'),
  plt,
  width = 40,
  height = 30,
  units = 'cm',
  bg = 'white',
  dpi = 600
)
