################################################################################
################################################################################

# DESCRIPTION:

# This script:
#   - filters plot level data
#   - tracks which plots/pixels are removed through filtering
#   - creates figures of plots/pixels removed through filtering

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
library(scales)
library(stringr)
library(lubridate)
library(dbplyr)

# 1.2 Parameters ---------------------------------------------------------------

in_version = 'v20240319'
out_version = 'v20240508'

n_plots_threshold = 3
n_pixels_threshold = 5
water_snow_threshold = 0.01 # Remove if more than 1% water/snow
rep_percent = 1 # Remove upper X% of data
rep_threshold = (100-rep_percent)/100
n_ccdc_obs_threshold = 30
doy_lwr = 'doy_startSnowfree'
doy_upr = 'doy_endSnowfree'
biomass_threshold = 15000
remove_built = TRUE

in_dir = 'output/04_predictor_extraction/'
out_dir = 'output/05_filtered/'
out_dir_secondary = 'output/06_model_ready_data/plot_level_mc/'
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

# 1.3 Read in data -------------------------------------------------------------

# Read in full synthesis dataset -- necessary for representativeness analysis
ds.all = fread(paste0(in_dir, 'arctic_tundra_all_biomass_synthesis_representativeness_predictors_', in_version, '.csv'))

for(ds_name in ds_name_list){
  
  # Read in current dataset
  ds = fread(paste0(in_dir, 'arctic_tundra_', ds_name, '_biomass_synthesis_representativeness_predictors_', in_version, '.csv'))
  
  # ==============================================================================
  # 2. TIDY DATA =================================================================
  # ==============================================================================
  
  # 2.1 Remove 'built' non-vegetated observations --------------------------------
  
  # Footprint of built areas often too small relative to Landsat pixel
  ds = ds[ds$site_id != 'built',]

  # 2.2 Create segment mismatch column -------------------------------------------
  
  # Calculate segment mismatch in days
  ds = ds %>%
    mutate(ccdc_segmentMismatchDays = case_when(ccdc_segmentType == 'current' ~ 0,
                                                ccdc_segmentType == 'previous' ~ as.numeric( difftime(parse_date_time(x = paste(year, as.character(harvest_doy)), orders = "yj"), date_decimal(ccdc_tEnd), units = 'days') ),
                                                ccdc_segmentType == 'next' ~ as.numeric( difftime(date_decimal(ccdc_tStart), parse_date_time(x = paste(year, as.character(harvest_doy)), orders = "yj"), units = 'days' ) ) ) ) %>%
    as.data.table()
  
  # 2.3 Create threshold value dataset -------------------------------------------
  
  ds.thresholds = data.frame(filter_name = c('ccdc_segment', 'water', 'snow', 'NBR', 'NDMI', 'NDVI', 'ccdc_obs', 'doy_lwr', 'doy_upr', 'tree_biomass', 'n_pixel', 'n_plot'), 
                             threshold_value = c(NA, as.character(water_snow_threshold), as.character(water_snow_threshold), as.character(rep_percent), as.character(rep_percent), 
                                                 as.character(rep_percent), as.character(n_ccdc_obs_threshold), 'doy_startSnowfree', 'doy_endSnowfree', as.character(biomass_threshold), 
                                                 as.character(n_pixels_threshold), as.character(n_plots_threshold)))
  
  # 2.4 Track original data and metadata -----------------------------------------
  
  # Count number of starting sites and plots, nonveg vs. veg

  # Assign variables -- start
  n_plot_start_nonveg = length(unique(ds[ds$dataset_id == 'nonveg',]$plot_code))
  n_plot_start_veg = length(unique(ds[ds$dataset_id != 'nonveg',]$plot_code))
  n_site_start_nonveg = length(unique(ds[ds$dataset_id == 'nonveg',]$site_code))
  n_site_start_veg = length(unique(ds[ds$dataset_id != 'nonveg',]$site_code))
  
  # Report
  print(paste0(ds_name, ' synthesis, original number of plots -- nonveg: ', n_plot_start_nonveg,
               ', veg: ', n_plot_start_veg))
  print(paste0(ds_name, ' synthesis, original number of sites -- nonveg: ', n_site_start_nonveg,
               ', veg: ', n_site_start_veg))
    
  # Retain copy of original datasets
  ds.orig = ds

  # ==============================================================================
  # 3. CHECK REPRESENTATIVENESS DATA =============================================
  # ==============================================================================
  
  # Coefficient of Variation not valid if both positive and negative data present
  # So, using standard deviation instead
  # First, need to ensure mean and standard deviation aren't too closely correlated
  # - https://stats.stackexchange.com/questions/56399/why-is-the-coefficient-of-variation-not-valid-when-using-data-with-positive-and
  # - https://personal.utdallas.edu/~herve/abdi-cv2010-pretty.pdf
  
  # 3.1 Check mean and standard deviation correlation ----------------------------
  
  nbr_cor = abs(cor(na.omit(ds.all$rep_NBR_mean), na.omit(ds.all$rep_NBR_stdDev)))
  ndmi_cor = abs(cor(na.omit(ds.all$rep_NDMI_mean), na.omit(ds.all$rep_NDMI_stdDev)))
  ndvi_cor = abs(cor(na.omit(ds.all$rep_NDVI_mean), na.omit(ds.all$rep_NDVI_stdDev)))
  
  print(paste0('NBR mean and standard deviation correlation: ', nbr_cor))
  print(paste0('NDMI mean and standard deviation correlation: ', ndmi_cor))
  print(paste0('NDVI mean and standard deviation correlation: ', ndvi_cor))
  if(nbr_cor >= 0.4){stop('NBR mean and standard deviation are moderately correlated (>= 0.4) - standard deviation might not be a good metric of representativeness')}
  if(ndmi_cor >= 0.4){stop('NDMI mean and standard deviation are moderately correlated (>= 0.4) - standard deviation might not be a good metric of representativeness')}
  if(ndvi_cor >= 0.4){stop('NDVI mean and standard deviation are moderately correlated (>= 0.4) - standard deviation might not be a good metric of representativeness')}
  
  nbr_cor = abs(cor(abs(na.omit(ds.all$rep_NBR_mean)), na.omit(ds.all$rep_NBR_stdDev)))
  ndmi_cor = abs(cor(abs(na.omit(ds.all$rep_NDMI_mean)), na.omit(ds.all$rep_NDMI_stdDev)))
  ndvi_cor = abs(cor(abs(na.omit(ds.all$rep_NDVI_mean)), na.omit(ds.all$rep_NDVI_stdDev)))
  
  print(paste0('NBR abs(mean) and standard deviation correlation: ', nbr_cor))
  print(paste0('NDMI abs(mean) and standard deviation correlation: ', ndmi_cor))
  print(paste0('NDVI abs(mean) and standard deviation correlation: ', ndvi_cor))
  if(nbr_cor >= 0.4){stop('NBR abs(mean) and standard deviation are moderately correlated (>= 0.4) - standard deviation might not be a good metric of representativeness')}
  if(ndmi_cor >= 0.4){stop('NDMI abs(mean) and standard deviation are moderately correlated (>= 0.4) - standard deviation might not be a good metric of representativeness')}
  if(ndvi_cor >= 0.4){stop('NDVI abs(mean) and standard deviation are moderately correlated (>= 0.4) - standard deviation might not be a good metric of representativeness')}
  
  # 3.2 Calculate representativeness quantiles ----------------------------
  
  # Perform on full synthesis dataset ('all') to keep filtering consistent among datasets
  nbr_std = quantile(na.omit(ds.all$rep_NBR_stdDev), rep_threshold)
  ndmi_std = quantile(na.omit(ds.all$rep_NDMI_stdDev), rep_threshold)
  ndvi_std = quantile(na.omit(ds.all$rep_NDVI_stdDev), rep_threshold)
  
  # ==============================================================================
  # 4. APPLY FILTERS =============================================================
  # ==============================================================================
  
  # 4.1 Identify high tree biomass sites -----------------------------------------
  
  # Identify sites that have any plot with tree biomass > threshold
  # We want to remove the entire site, otherwise we are just eliminating plots with high tree biomass, but later averaging over the remaining plots
  # This will give a poor representation of woody biomass at the site
  high_biomass_sites = ds.all[ds.all$pft =='tree' & ds.all$biomass_density_gm2 > biomass_threshold,]$site_code
  
  # 4.2 Loop datasets and filter -------------------------------------------------
  
  # CCDC segment mismatch
  ds = ds[ds$ccdc_segmentMismatchDays <= 120,]
  
  # Water and snow 
  ds = ds[ds$rep_snow_fraction <= water_snow_threshold | ds$dataset_id == 'nonveg',]
  ds = ds[ds$rep_water_fraction <= water_snow_threshold | ds$dataset_id == 'nonveg',]
  
  # CCDC observations 
  ds = ds[ds$ccdc_numObs >= n_ccdc_obs_threshold,]
  
  # Valid spectral band data
  ds = ds %>%
    filter_at(vars(matches("red|green|blue|NIR|SWIR1|SWIR2") & !matches("ccdc|coefs|change|Range|rep|texture|removed")),
              all_vars(.>=0 & .<=10000))
  
  # Valid spectral index data
  ds = ds %>%
    filter_at(vars(matches("EVI2b|NBR|NDMI|NDVI|NDWI") & !matches("ccdc|coefs|change|Range|rep|texture|removed")),
              any_vars(.>=-10000 & .<=10000))
  
  # No representativeness data 
  ds = ds[!is.na(ds$rep_NBR_mean),]
  
  # NBR representativeness
  ds = ds[ds$rep_NBR_stdDev <= nbr_std,]
  
  # NDMI representativeness
  ds = ds[ds$rep_NDMI_stdDev <= ndmi_std,]
  
  # NDVI representativeness
  ds = ds[ds$rep_NDVI_stdDev <= ndvi_std,]
  
  # Day-of-year
  ds = ds[(ds$harvest_doy >= ds[[doy_lwr]] & ds$harvest_doy <= ds[[doy_upr]]) | ds$dataset_id == 'nonveg',]
  
  # High tree biomass
  ds = ds[!(ds$site_code %in% high_biomass_sites),]
  
  # Calculate number of plots and pixels per site
  if(ds_name == 'all'){
    ds = ds %>%
      group_by(site_code, pft) %>%
      mutate(
        n_plots = n_distinct(plot_code),
        n_pixels = n_distinct(pixel_num)
      ) %>% 
      as.data.table()
  }else{
    ds = ds %>%
      group_by(site_code) %>%
      mutate(
        n_plots = n_distinct(plot_code),
        n_pixels = n_distinct(pixel_num)
      ) %>% 
      as.data.table()
  }
  # Save version after applying all filters but n_plot/n_pixel for use in tracking removed data
  ds.filtered = ds
  
  # Remove any site level observations with fewer than 5 pixels per site (max is 9)
  ds = ds[ds$n_pixels >= n_pixels_threshold | ds$coord_type == 'plot' | ds$dataset_id == 'nonveg',] 
  
  # Remove any observations with fewer than 3 plots per site
  # NOTE: ds25 site biomass estimates were calculated from 5-10 subplots, but provided to us in aggregated form thus they only have 1 plot per site -- exempt these from the >=3 plots per site filter
  ds = ds[ds$n_plots >= n_plots_threshold | ds$dataset_id == 'ds25' | ds$dataset_id == 'nonveg',] 
  
  # ==============================================================================
  # 5. SAVE FILTERED DATA ========================================================
  # ==============================================================================
  
  print(paste0('Final number of plots for all dataset -- nonveg: ', length(unique(ds[ds$dataset_id == 'nonveg',]$plot_code)), ', veg: ', length(unique(ds[ds$dataset_id != 'nonveg',]$plot_code))))
  print(paste0('Final number of sites for all dataset -- nonveg: ', length(unique(ds[ds$dataset_id == 'nonveg',]$site_code)), ', veg: ', length(unique(ds[ds$dataset_id != 'nonveg',]$site_code))))

  fwrite(ds, paste0(out_dir, 'ds_', ds_name, '_plots_filtered_', out_version, '.csv'))
  fwrite(ds, paste0(out_dir_secondary, 'ds_', ds_name, '_plots_filtered_', out_version, '.csv'))
  
  fwrite(ds.thresholds, paste0(out_dir, 'ds_thresholds_', out_version, '.csv'))
  
  # ==============================================================================
  # 6. TRACK REMOVED DATA ========================================================
  # ==============================================================================
  
  # 6.1 Track removed data -----------------------------------------------------
  
  # Get data
  ds.removed = ds.orig
  
  # Track removed data
  ds.removed$ccdc_segment_removed = with(ds.orig, ifelse(ccdc_segmentMismatchDays <= 120, FALSE, TRUE))
  ds.removed$water_removed = with(ds.orig, ifelse(rep_water_fraction <= water_snow_threshold | dataset_id == 'nonveg', FALSE, TRUE))
  ds.removed$snow_removed = with(ds.orig, ifelse(rep_snow_fraction <= water_snow_threshold | dataset_id == 'nonveg', FALSE, TRUE))
  ds.removed$ccdc_obs_removed = with(ds.orig, ifelse(ccdc_numObs >= n_ccdc_obs_threshold, FALSE, TRUE))
  ds.removed$no_rep_removed = with(ds.orig, ifelse(!is.na(rep_NBR_mean), FALSE, TRUE))
  ds.removed$nbr_removed = with(ds.orig, ifelse(rep_NBR_stdDev <= nbr_std, FALSE, TRUE))
  ds.removed$ndmi_removed = with(ds.orig, ifelse(rep_NDMI_stdDev <= ndmi_std, FALSE, TRUE))
  ds.removed$ndvi_removed = with(ds.orig, ifelse(rep_NDVI_stdDev <= ndvi_std, FALSE, TRUE))
  ds.removed$doy_removed = with(ds.orig, ifelse((harvest_doy >= ds.orig[[doy_lwr]] & harvest_doy <= ds.orig[[doy_upr]]) | dataset_id == 'nonveg', FALSE, TRUE))
  ds.removed$high_bmass_removed = with(ds.orig, ifelse(!(site_code %in% high_biomass_sites), FALSE, TRUE))
  
  # Track removed data - adjust for spectral validity
  ds.removed.keep = ds.removed %>%
    filter_at(vars(matches("red|green|blue|NIR|SWIR1|SWIR2") & !matches("ccdc|coefs|change|Range|rep|texture|removed")),
              all_vars(.>=0 & .<=10000))
  ds.removed.remove = ds.removed %>%
    filter_at(vars(matches("red|green|blue|NIR|SWIR1|SWIR2") & !matches("ccdc|coefs|change|Range|rep|texture|removed")),
              any_vars(.<0 | .>10000))
  ds.removed.keep$spectral_removed = FALSE
  ds.removed.remove$spectral_removed = TRUE
  ds.removed = rbind(ds.removed.keep, ds.removed.remove)
  
  # Track removed data - adjust for index validity
  ds.removed.keep = ds.removed %>%
    filter_at(vars(matches("EVI2b|NBR|NDMI|NDVI|NDWI") & !matches("ccdc|coefs|change|Range|rep|texture|removed")),
              all_vars(.>= -10000 & .<=10000))
  ds.removed.remove = ds.removed %>%
    filter_at(vars(matches("EVI2b|NBR|NDMI|NDVI|NDWI") & !matches("ccdc|coefs|change|Range|rep|texture|removed")),
              any_vars(.< -10000 | .>10000))
  ds.removed.keep$index_removed = FALSE
  ds.removed.remove$index_removed = TRUE
  ds.removed = rbind(ds.removed.keep, ds.removed.remove)
  
  # For n_plots and n_pixels, we want to perform filter after all other filtering steps have been completed
  # ds.filtered has fewer rows than ds.removed so we need peform a series of steps
  
  # Assign yes/no on ds.filtered
  ds.filtered$n_pixels_removed = with(ds.filtered, ifelse(n_pixels >= n_pixels_threshold | coord_type == 'plot' | dataset_id == 'nonveg', FALSE, TRUE))
  ds.filtered$n_plots_removed = with(ds.filtered, ifelse(n_plots >= n_plots_threshold | dataset_id == 'ds25' | dataset_id == 'nonveg', FALSE, TRUE))
  
  # Join with ds.removed
  if(ds_name == 'all'){
    ds.removed = dplyr::left_join(ds.removed, ds.filtered %>% select(c('site_code', 'plot_code', 'pixel_num', 'n_pixels', 'n_plots', 'n_pixels_removed', 'n_plots_removed', 'pft')), by = c('site_code', 'plot_code', 'pixel_num', 'pft')) %>% as.data.table()
  }else{
    ds.removed = dplyr::left_join(ds.removed, ds.filtered %>% select(c('site_code', 'plot_code', 'pixel_num', 'n_pixels', 'n_plots', 'n_pixels_removed', 'n_plots_removed')), by = c('site_code', 'plot_code', 'pixel_num')) %>% as.data.table()
  }  
  
  # Fill NAs with FALSE
  ds.removed$n_pixels_removed[is.na(ds.removed$n_pixels_removed)] = FALSE
  ds.removed$n_plots_removed[is.na(ds.removed$n_plots_removed)] = FALSE
  
  # Format original and filtered data
  ds.orig$veg_nonveg = 'Vegetated'
  ds.orig$veg_nonveg[ds.orig$dataset_id == 'nonveg'] = 'Non-vegetated'
  ds.orig$veg_nonveg = factor(ds.orig$veg_nonveg, levels=c('Vegetated', 'Non-vegetated'))
  ds.filtered$veg_nonveg = 'Vegetated'
  ds.filtered$veg_nonveg[ds.filtered$dataset_id == 'nonveg'] = 'Non-vegetated'
  ds.filtered$veg_nonveg = factor(ds.filtered$veg_nonveg, levels=c('Vegetated', 'Non-vegetated'))
  
  # 6.2 Plot all data ----------------------------------------------------------
  
  ds.orig = data.frame(ds.orig)
  ds.filtered = data.frame(ds.filtered)
  
  # CCDC segment
  ccdc_segment = ggplot(arrange(ds.orig, ccdc_segmentMismatchDays), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = ccdc_segmentMismatchDays, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[ds.orig$ccdc_segmentMismatchDays > 120,], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'CCDC\nSegment\nMismatch\nDays', title = 'Remove if > 120 day mismatch', shape = '')+
    theme_minimal(base_size = 30)
  
  # Water
  water = ggplot(arrange(ds.orig, rep_water_fraction), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = rep_water_fraction, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[ds.orig$rep_water_fraction > water_snow_threshold & ds.orig$dataset_id != 'nonveg',], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'Water Fraction', title = paste0('Remove if > ', water_snow_threshold), shape = '')+
    theme_minimal(base_size = 30)
  
  # Snow
  snow = ggplot(arrange(ds.orig, rep_snow_fraction), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = rep_snow_fraction, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[ds.orig$rep_snow_fraction > water_snow_threshold & ds.orig$dataset_id != 'nonveg',], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'Snow Fraction', title = paste0('Remove if > ', water_snow_threshold), shape = '')+
    theme_minimal(base_size = 30)
  
  # CCDC observations
  ccdc_obs = ggplot(arrange(ds.orig, ccdc_numObs), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = ccdc_numObs, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[ds.orig$ccdc_numObs < n_ccdc_obs_threshold,], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = '# CCDC\nObservations', title = paste0('Remove if < ', n_ccdc_obs_threshold), shape = '')+
    theme_minimal(base_size = 30)
  
  # Valid spectral data
  spectral = ggplot(arrange(ds.orig, spectral_blue_lateSummer), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = spectral_blue_lateSummer, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig %>% filter_at(vars(matches("EVI2b|NBR|NDMI|NDVI|NDWI") & !matches("ccdc|coefs|change|Range|rep|texture|removed")), any_vars(.< -10000 | .>10000)), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'blue\nlate summer\n(ex. band)', title = 'Remove if any spectral data < 0 or > 10,000', shape = '')+
    theme_minimal(base_size = 30)
  
  # Valid index data
  index = ggplot(arrange(ds.orig, spectral_NDWI_lateSummer), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = spectral_NDWI_lateSummer, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig %>% filter_at(vars(matches("EVI2b|NBR|NDMI|NDVI|NDWI") & !matches("ccdc|coefs|change|Range|rep|texture|removed")), any_vars(.< -10000 | .>10000)), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'NDWI\nlate summer\n(ex. band)', title = 'Remove if any index data < -10,000 or > 10,000', shape = '')+
    theme_minimal(base_size = 30)
  
  # No representativeness data
  no_rep = ggplot(ds.orig, aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[is.na(ds.orig$rep_NBR_mean),], col = 'red', size = 2)+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', title = 'Remove if no representativeness data', shape = '')+
    theme_minimal(base_size = 30)
  
  # NBR
  nbr = ggplot(arrange(ds.orig, rep_NBR_stdDev), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = rep_NBR_stdDev, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = na.omit(ds.orig[ds.orig$rep_NBR_stdDev > nbr_std,]), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'NBR\nStdDev', title = paste0('Remove if >= ', round(nbr_std, 0)), shape = '')+
    theme_minimal(base_size = 30)
  
  # NDMI
  ndmi = ggplot(arrange(ds.orig, rep_NDMI_stdDev), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = rep_NDMI_stdDev, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = na.omit(ds.orig[ds.orig$rep_NDMI_stdDev > ndmi_std,]), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'NDMI\nStdDev', title = paste0('Remove if >= ', round(ndmi_std, 0)), shape = '')+
    theme_minimal(base_size = 30)
  
  # NDVI
  ndvi = ggplot(arrange(ds.orig, rep_NDVI_stdDev), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = rep_NDVI_stdDev, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = na.omit(ds.orig[ds.orig$rep_NDVI_stdDev > ndvi_std,]), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'NDVI\nStdDev', title = paste0('Remove if >= ', round(ndvi_std, 0)), shape = '')+
    theme_minimal(base_size = 30)
  
  # Harvest day of year
  doy = ggplot(ds.orig, aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = harvest_doy, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[(ds.orig$harvest_doy < ds.orig[[doy_lwr]] | ds.orig$harvest_doy > ds.orig[[doy_upr]]),], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'Harvest\nDOY', title = paste0('Remove if < ', gsub('_doy', '', doy_lwr), ' or > ', gsub('_doy', '', doy_upr)), shape = '')+
    theme_minimal(base_size = 30)
  
  # High tree biomass
  woody_biomass = ggplot(ds.orig, aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[(ds.orig$site_code %in% high_biomass_sites),], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', title = paste0('Remove if tree biomass > ', biomass_threshold), shape = '')+
    theme_minimal(base_size = 30)
  
  # Number of pixels
  n_pixels = ggplot(arrange(ds.filtered, desc(n_pixels)), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = n_pixels, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = na.omit(ds.filtered[ds.filtered$n_pixels < n_pixels_threshold & ds.filtered$coord_type != 'plot' & ds.filtered$dataset_id != 'nonveg',]), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = '# of Pixels', title = paste0('Remove if < ', round(n_pixels_threshold, 0)), shape = '')+
    theme_minimal(base_size = 30)
  
  # Number of plots
  n_plots = ggplot(arrange(ds.filtered, desc(n_plots)), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = n_plots, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = na.omit(ds.filtered[ds.filtered$n_plots < n_plots_threshold & ds.filtered$dataset_id != 'ds25' & ds.filtered$dataset_id != 'nonveg',]), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = '# of Plots', title = paste0('Remove if < ', round(n_plots_threshold, 0)), shape = '')+
    theme_minimal(base_size = 30)
  
  # 6.2 Plot vegetated data only -----------------------------------------------
  
  ds.orig = ds.orig[ds.orig$dataset_id != 'nonveg',]
  ds.filtered = ds.filtered[ds.filtered$dataset_id != 'nonveg',]
  
  # CCDC segment
  ccdc_segment_veg = ggplot(arrange(ds.orig, ccdc_segmentMismatchDays), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = ccdc_segmentMismatchDays, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[ds.orig$ccdc_segmentMismatchDays > 120,], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'CCDC\nSegment\nMismatch\nDays', title = 'Remove if > 120 day mismatch', shape = '')+
    theme_minimal(base_size = 30)
  
  # Water
  water_veg = ggplot(arrange(ds.orig, rep_water_fraction), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = rep_water_fraction, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[ds.orig$rep_water_fraction > water_snow_threshold & ds.orig$dataset_id != 'nonveg',], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'Water Fraction', title = paste0('Remove if > ', water_snow_threshold), shape = '')+
    theme_minimal(base_size = 30)
  
  # Snow
  snow_veg = ggplot(arrange(ds.orig, rep_snow_fraction), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = rep_snow_fraction, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[ds.orig$rep_snow_fraction > water_snow_threshold & ds.orig$dataset_id != 'nonveg',], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'Snow Fraction', title = paste0('Remove if > ', water_snow_threshold), shape = '')+
    theme_minimal(base_size = 30)
  
  # CCDC observations
  ccdc_obs_veg = ggplot(arrange(ds.orig, ccdc_numObs), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = ccdc_numObs, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[ds.orig$ccdc_numObs < n_ccdc_obs_threshold,], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = '# CCDC\nObservations', title = paste0('Remove if < ', n_ccdc_obs_threshold), shape = '')+
    theme_minimal(base_size = 30)
  
  # Valid spectral data
  spectral_veg = ggplot(arrange(ds.orig, spectral_blue_lateSummer), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = spectral_blue_lateSummer, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig %>% filter_at(vars(matches("EVI2b|NBR|NDMI|NDVI|NDWI") & !matches("ccdc|coefs|change|Range|rep|texture|removed")), any_vars(.< -10000 | .>10000)), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'blue\nlate summer\n(ex. band)', title = 'Remove if any spectral data < 0 or > 10,000', shape = '')+
    theme_minimal(base_size = 30)
  
  # Valid index data
  index_veg = ggplot(arrange(ds.orig, spectral_NDWI_lateSummer), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = spectral_NDWI_lateSummer, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig %>% filter_at(vars(matches("EVI2b|NBR|NDMI|NDVI|NDWI") & !matches("ccdc|coefs|change|Range|rep|texture|removed")), any_vars(.< -10000 | .>10000)), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'NDWI\nlate summer\n(ex. band)', title = 'Remove if any index data < -10,000 or > 10,000', shape = '')+
    theme_minimal(base_size = 30)
  
  # No representativeness data
  no_rep_veg = ggplot(ds.orig, aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[is.na(ds.orig$rep_NBR_mean),], col = 'red', size = 2)+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', title = 'Remove if no representativeness data', shape = '')+
    theme_minimal(base_size = 30)
  
  # NBR
  nbr_veg = ggplot(arrange(ds.orig, rep_NBR_stdDev), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = rep_NBR_stdDev, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = na.omit(ds.orig[ds.orig$rep_NBR_stdDev > nbr_std,]), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'NBR\nStdDev', title = paste0('Remove if >= ', round(nbr_std, 0)), shape = '')+
    theme_minimal(base_size = 30)
  
  # NDMI
  ndmi_veg = ggplot(arrange(ds.orig, rep_NDMI_stdDev), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = rep_NDMI_stdDev, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = na.omit(ds.orig[ds.orig$rep_NDMI_stdDev > ndmi_std,]), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'NDMI\nStdDev', title = paste0('Remove if >= ', round(ndmi_std, 0)), shape = '')+
    theme_minimal(base_size = 30)
  
  # NDVI
  ndvi_veg = ggplot(arrange(ds.orig, rep_NDVI_stdDev), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = rep_NDVI_stdDev, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = na.omit(ds.orig[ds.orig$rep_NDVI_stdDev > ndvi_std,]), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'NDVI\nStdDev', title = paste0('Remove if >= ', round(ndvi_std, 0)), shape = '')+
    theme_minimal(base_size = 30)
  
  # Harvest day of year
  doy_veg = ggplot(ds.orig, aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = harvest_doy, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[(ds.orig$harvest_doy < ds.orig[[doy_lwr]] | ds.orig$harvest_doy > ds.orig[[doy_upr]]),], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = 'Harvest\nDOY', title = paste0('Remove if < ', gsub('_doy', '', doy_lwr), ' or > ', gsub('_doy', '', doy_upr)), shape = '')+
    theme_minimal(base_size = 30)
  
  # High woody biomass
  woody_biomass_veg = ggplot(ds.orig, aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = ds.orig[(ds.orig$site_code %in% high_biomass_sites),], col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', title = paste0('Remove if tree biomass > ', biomass_threshold), shape = '')+
    theme_minimal(base_size = 30)
  
  # Number of pixels
  n_pixels_veg = ggplot(arrange(ds.filtered, desc(n_pixels)), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = n_pixels, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = na.omit(ds.filtered[ds.filtered$n_pixels < n_pixels_threshold & ds.filtered$coord_type != 'plot' & ds.filtered$dataset_id != 'nonveg',]), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = '# of Pixels', title = paste0('Remove if < ', round(n_pixels_threshold, 0)), shape = '')+
    theme_minimal(base_size = 30)
  
  # Number of plots
  n_plots_veg = ggplot(arrange(ds.filtered, desc(n_plots)), aes(x = spectral_NDVI_peakSummer, y = biomass_density_gm2, col = n_plots, shape = veg_nonveg))+
    geom_point(size = 4)+
    geom_point(data = na.omit(ds.filtered[ds.filtered$n_plots < n_plots_threshold & ds.filtered$dataset_id != 'ds25' & ds.filtered$dataset_id != 'nonveg',]), col = 'red', size = 2)+
    scale_color_continuous(type = "viridis")+
    scale_shape_manual(labels = c("Vegetated","Non-vegetated"),
                       values = c(16, 17),
                       guide = guide_legend(override.aes =list(col = "black")))+
    labs(y = paste0(str_to_title(ds_name), ' Plant Biomass (g/m2)'), x = 'Peak Summer NDVI', col = '# of Plots', title = paste0('Remove if < ', round(n_plots_threshold, 0)), shape = '')+
    theme_minimal(base_size = 30)
  
  # 6.3 Save plots -------------------------------------------------------------
  
  plts = list(ccdc_segment,
              ccdc_segment_veg,
              water, 
              water_veg,
              snow, 
              snow_veg,
              ccdc_obs, 
              ccdc_obs_veg,
              spectral,
              spectral_veg,
              index,
              index_veg,
              no_rep, 
              no_rep_veg,
              nbr, 
              nbr_veg,
              ndmi, 
              ndmi_veg, 
              ndvi, 
              ndvi_veg, 
              doy, 
              doy_veg, 
              woody_biomass, 
              woody_biomass_veg,
              n_pixels,
              n_pixels_veg,
              n_plots,
              n_plots_veg)
  
  plt_names = list(paste0('ccdc_segment_', ds_name),
                   paste0('ccdc_segment_veg_', ds_name),
                   paste0('water_', ds_name),
                   paste0('water_veg_', ds_name),
                   paste0('snow_', ds_name),
                   paste0('snow_veg_', ds_name),
                   paste0('ccdc_obs_', ds_name),
                   paste0('ccdc_obs_veg_', ds_name),
                   paste0('spectral_', ds_name),
                   paste0('spectral_veg_', ds_name),
                   paste0('index_', ds_name),
                   paste0('index_veg_', ds_name),
                   paste0('no_rep_', ds_name),
                   paste0('no_rep_veg_', ds_name),
                   paste0('nbr_', ds_name),
                   paste0('nbr_veg_', ds_name),
                   paste0('ndmi_', ds_name),
                   paste0('ndmi_veg_', ds_name),
                   paste0('ndvi_', ds_name),
                   paste0('ndvi_veg_', ds_name),
                   paste0('doy_', ds_name),
                   paste0('doy_veg_', ds_name),
                   paste0('woody_biomass_', ds_name),
                   paste0('woody_biomass_veg_', ds_name),
                   paste0('n_pixels_', ds_name),
                   paste0('n_pixels_veg_', ds_name),
                   paste0('n_plots_', ds_name),
                   paste0('n_plots_veg_', ds_name))
  
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
  # 7. SAVE REMOVED DATA =========================================================
  # ==============================================================================
  
  fwrite(ds.removed, paste0(out_dir, 'ds_', ds_name, '_plots_removed_', out_version, '.csv'), row.names = FALSE)

}