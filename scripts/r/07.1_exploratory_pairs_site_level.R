################################################################################
################################################################################

# DESCRIPTION:

# This script produces pairs plots of response (biomass) and predictors

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

# NOTES:

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

rm(list=ls())
require(data.table)
require(tidyverse)
require(scales)

# 1.2 Parameters ---------------------------------------------------------------

in_dir = 'output/06_model_ready_data/'
version = 'v20240508'
response_type = 'all_data' # Choose 'all_data' or 'presence_only'
out_fig_dir = 'output/07_exploratory_analysis/pairs/'

# Check if response type subdirectory exists, if not create it
if(file.exists(out_fig_dir)){
  print('Output version directory exists')
  cat('\n')
} else {
  dir.create(out_fig_dir)
  print('Output version directory created')
  cat('\n')
}

# 1.3 Read in data -------------------------------------------------------------

ds.woody = fread(paste0(in_dir, 'ds_woody_', version, '.csv'))
ds.total = fread(paste0(in_dir, 'ds_total_', version, '.csv'))

# ==============================================================================
# 2. PLOT ======================================================================
# ==============================================================================

# 2.1 Set categorical variables ------------------------------------------------

ds.woody$zone_3 = as.factor(ds.woody$zone_3)
ds.woody$world_terrestrial_ecosystems = as.factor(ds.woody$world_terrestrial_ecosystems)
ds.woody$ecoregion = as.factor(ds.woody$ecoregion)

ds.total$zone_3 = as.factor(ds.total$zone_3)
ds.total$world_terrestrial_ecosystems = as.factor(ds.total$world_terrestrial_ecosystems)
ds.total$ecoregion = as.factor(ds.total$ecoregion)

# 2.2 Tidy ------------------------------------------------

if(response_type == 'presence_only'){
  ds.woody = ds.woody[ds.woody$plant_biomass_density_gm2 > 0,]
  ds.total = ds.total[ds.total$plant_biomass_density_gm2 > 0,]
}

predictor_names = grep("permafrost|spectral_|texture_|topo_|trees_|zone_3|ecoregion|world_terrestrial_ecosystems|latitude|longitude", names(ds.total), value = TRUE)

# 2.3 Pairs plots ------------------------------------------------

i = 1
for(predictor_name in predictor_names){
  
  print(paste0('Woody biomass, exporting figure ', i, ' of ', length(predictor_names)))
  i = i + 1
        
  plt = ggplot(ds.woody, aes(x = !!sym(predictor_name), y = biomass_density_gm2, color = zone_3))+
    geom_point(size = 3)+
    labs(x = predictor_name, y = 'Woody Biomass (g/m2)', color = 'Zone')+
    theme_minimal(base_size = 30)
  plt
  
  ggsave(
    paste0(out_fig_dir, 'woody_biomass_', predictor_name, '.png'),
    plt,
    width = 40,
    height = 30,
    units = 'cm',
    bg = 'white',
    dpi = 600
  )
  
}

i = 1
for(predictor_name in predictor_names){
  
  print(paste0('Total biomass, exporting figure ', i, ' of ', length(predictor_names)))
  i = i + 1
  
  plt = ggplot(ds.total, aes(x = !!sym(predictor_name), y = biomass_density_gm2, color = zone_3))+
    geom_point(size = 3)+
    labs(x = predictor_name, y = 'Total Plant Biomass (g/m2)', color = 'Zone')+
    theme_minimal(base_size = 30)
  plt
  
  ggsave(
    paste0(out_fig_dir, 'total_biomass_', predictor_name, '.png'),
    plt,
    width = 40,
    height = 30,
    units = 'cm',
    bg = 'white',
    dpi = 600
  )
  
}

