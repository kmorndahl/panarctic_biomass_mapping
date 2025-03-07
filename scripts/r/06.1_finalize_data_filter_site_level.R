################################################################################
################################################################################

# DESCRIPTION:

# This script:
#   - filters site level data
#   - tracks which sites are removed through filtering
#   - creates figures of sites removed through filtering

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
require(sf)

# 1.2 Parameters ---------------------------------------------------------------

version = 'v20240508'

plot_area_threshold = 0.05
ccdc_rmse_threshold = 2500

filter_df = data.frame(filter_type = character(), filter_threshold = integer(), dataset_type = character(), land_cover = character(), n_removed = integer())

in_dir = 'output/05_filtered/'
out_dir = 'output/06_model_ready_data/'
out_fig_dir = 'output/07_exploratory_analysis/filtering/'

ds_name_list = c('total', 'woody')

# Check if response type subdirectory exists, if not create it
if(file.exists(out_fig_dir)){
  print('Output version directory exists')
  cat('\n')
} else {
  dir.create(out_fig_dir)
  print('Output version directory created')
  cat('\n')
}

for(ds_name in ds_name_list){
  
  # 1.3 Read in data -------------------------------------------------------------
  
  ds = read.csv(paste0(in_dir, 'ds_', ds_name, '_site_level_', version, '.csv'))
  
  # ==============================================================================
  # 2. APPLY FILTERS =============================================================
  # ==============================================================================
  
  # 2.1 Tidy -------------------------------------------------------------------
  
  # Calculate average CCDC RMSE across all bands/indices
  ds = ds %>% mutate(ds, ccdc_rmse_avg = rowMeans(select(ds, ends_with("_rmse")), na.rm = TRUE)) 
  
  # Report starting number of sites
  print(paste0(str_to_title(ds_name), ' synthesis, starting number of sites -- nonveg: ', length(unique(ds[ds$dataset_id == 'nonveg',]$site_code)), ', veg: ', length(unique(ds[ds$dataset_id != 'nonveg',]$site_code))))
  
  # Save copy of unfiltered data
  ds.orig = ds
  
  # 2.2 Apply filters ---------------------------------------------------------
  
  # Plot area
  ds = ds[ds$plot_area_m2_mean >= plot_area_threshold,]
  
  # CCDC RMSE
  ds = ds[ds$ccdc_rmse_avg < ccdc_rmse_threshold,]
  
  # 2.3 Track removed data -----------------------------------------------------
  
  ds.removed = ds.orig
  
  ds.removed$plot_area_removed = with(ds.orig, ifelse(plot_area_m2_mean >= plot_area_threshold, FALSE, TRUE))
  ds.removed$ccdc_rmse_removed = with(ds.orig, ifelse(ccdc_rmse_avg < ccdc_rmse_threshold, FALSE, TRUE))
  
  # 2.4 Plot -----------------------------------------------------
  
  # All data -----
  
  ds.orig$veg_nonveg = 'Vegetated'
  ds.orig$veg_nonveg[ds.orig$dataset_id == 'nonveg'] = 'Non-vegetated'
  ds.orig$veg_nonveg = factor(ds.orig$veg_nonveg, levels=c('Vegetated', 'Non-vegetated'))
  
  # Plot size
  plot_size = ggplot(data.frame(arrange(ds.orig, desc(plot_area_m2_mean))), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = plot_area_m2_mean, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[ds.orig$plot_area_m2_mean < plot_area_threshold,], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis", limits = c(0,0.25), oob = scales::squish)+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'Plot Area (m2)', title = paste0('Remove if < ', plot_area_threshold))+
    theme_minimal(base_size = 30)
  
  # CCDC RMSE
  ccdc_rmse = ggplot(data.frame(arrange(ds.orig, ccdc_rmse_avg)), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = ccdc_rmse_avg, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[ds.orig$ccdc_rmse_avg >= ccdc_rmse_threshold,], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis", limits = c(0, max(ds.orig$ccdc_rmse_avg)), oob = scales::squish)+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'Avg. CCDC RMSE', title = paste0('Remove if any >= ', ccdc_rmse_threshold))+
    theme_minimal(base_size = 30)
  
  # Vegetated only -----
  
  ds.orig = ds.orig[ds.orig$dataset_id != 'nonveg',]
  
  # Plot size
  plot_size_veg = ggplot(data.frame(arrange(ds.orig, desc(plot_area_m2_mean))), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = plot_area_m2_mean))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[ds.orig$plot_area_m2_mean < plot_area_threshold,], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis", limits = c(0,0.25), oob = scales::squish)+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'Plot Area (m2)', title = paste0('Remove if < ', plot_area_threshold))+
    theme_minimal(base_size = 30)
  
  # CCDC RMSE
  ccdc_rmse_veg = ggplot(data.frame(arrange(ds.orig, ccdc_rmse_avg)), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = ccdc_rmse_avg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[ds.orig$ccdc_rmse_avg >= ccdc_rmse_threshold,], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis", limits = c(0, max(ds.orig$ccdc_rmse_avg)), oob = scales::squish)+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'Avg. CCDC RMSE', title = paste0('Remove if any >= ', ccdc_rmse_threshold))+
    theme_minimal(base_size = 30)
  
  # Name and export -----
  
  plts = list(plot_size,
              plot_size_veg,
              ccdc_rmse, 
              ccdc_rmse_veg)
  
  plt_names = list(paste0('plot_size_', ds_name),
                   paste0('plot_size_veg_', ds_name),
                   paste0('ccdc_rmse_', ds_name),
                   paste0('ccdc_rmse_veg_', ds_name))
  
  plts = setNames(plts, plt_names)
  
  for(i in 1:length(plts)){
    
    out_name = names(plts[i])
    
    ggsave(
      paste0(out_fig_dir, out_name, '.png'),
      plts[[i]],
      width = 40,
      height = 30,
      units = 'cm',
      bg = 'white',
      dpi = 600
    )
    
  }
  
  # ==============================================================================
  # 3. SAVE ======================================================================
  # ==============================================================================
  
  # 3.1 Filtered data ------------------------------------------------------------
  
  # Report
  print(paste0('Final number of sites for dataset -- nonveg: ', length(unique(ds[ds$dataset_id == 'nonveg',]$site_code)), ', veg: ', length(unique(ds[ds$dataset_id != 'nonveg',]$site_code))))
  
  # Save
  write.csv(ds, paste0(out_dir, 'ds_', ds_name, '_', version, '.csv'), row.names = FALSE)
  
  # 3.2 Removed data -------------------------------------------------------------
  
  # Save
  write.csv(ds.removed, paste0(in_dir, 'ds_', ds_name, '_sites_removed_', version, '.csv'), row.names = FALSE)
  
}