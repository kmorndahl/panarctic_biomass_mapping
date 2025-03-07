################################################################################
################################################################################

# DESCRIPTION:

# This script aggregates plot level data to the site level

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

# NOTES:
# - https://stackoverflow.com/questions/63382926/mutate-with-across-and-names-glue-cannot-interpolate-functions-into-strings

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

rm(list=ls())
require(data.table)
require(tidyverse)

# 1.2 Parameters ---------------------------------------------------------------

in_dir = 'output/05_filtered/'
out_dir = 'output/05_filtered/'

version = 'v20240508'

ds_name_list = c('total', 'woody')

# 1.3 Functions -------------------------------------------------------------

find_mode <- function(x) {
  u <- unique(na.omit(x))
  tab <- tabulate(match(x, u))
  u[tab == max(tab)][1] # If there are two modes, grab the first one
}

for(ds_name in ds_name_list){
  
  # 1.4 Read in data -------------------------------------------------------------
  
  ds = fread(paste0(in_dir, 'ds_', ds_name, '_plots_filtered_', version, '.csv'))
  
  # ==============================================================================
  # 2. AGGREGATE =================================================================
  # ==============================================================================
  
  # 2.1 Prepare dataset for aggregation -----------------------------------------
  
  # Specify columns to remove
  remove_cols = grep("biomass_dry_weight_g|rep_|pixel_num|distance|n_plots|n_pixels|ccdc_segment|ccdc_nSegments|ccdc_next|ccdc_previous|ccdc_t", names(ds), value = TRUE)
  
  # 2.2 Aggregate ----------------------------------------------------------------
  
  # Remove unneeded columns
  ds = ds %>% select(-any_of(remove_cols)) %>% as.data.table()
  
  # Get columns to average
  summarise_cols = grep("ccdc_|coefs_|doy_|spectral_|texture_|climate_|topo_|soil_|trees_cover|permafrost|biomass_density_gm2|plot_area_m2_mean|latitude|longitude", names(ds), value = TRUE)
  
  # Convert columns to numeric if they are to be averaged
  # There will be some NAs introduced for: ccdc_NextSegmentStart, ccdc_NextSegmentEnd, ccdc_previousSegmentStart, ccdc_previousSegmentEnd
  # That is okay, this is where these are 'no segment'
  ds = ds %>% mutate_at(summarise_cols, as.numeric) %>% as.data.table()
  
  # A few NAs for soil and world_terrestrial_ecosystems
  # These are nonveg points in ice or water, we can safely remove them
  drop_na_cols = grep("soil_|world_terrestrial_ecosystems", names(ds), value = TRUE)
  ds = ds %>% drop_na(any_of(drop_na_cols)) %>% as.data.table()
  
  # Report starting number of sites
  n_start = length(unique(ds$site_code))
  print(paste0(str_to_title(ds_name), ' synthesis, starting number of sites: ', n_start))
  
  # Aggregate data
  if(ds_name == 'all'){
    ds.agg = ds %>%
      group_by(site_code, pft) %>% # Also group by pft for 'all' dataset
      summarise(
        across(all_of(summarise_cols), mean, na.rm = TRUE), # For numeric columns, get mean
        across(-summarise_cols, find_mode) # For non-numeric columns, get mode
      ) %>%
      as.data.table()
  }else{
    ds.agg = ds %>%
      group_by(site_code) %>%
      summarise(
        across(all_of(summarise_cols), mean, na.rm = TRUE), # For numeric columns, get mean
        across(-summarise_cols, find_mode) # For non-numeric columns, get mode
      ) %>%
      as.data.table()
  }
  
  # ==============================================================================
  # 3. CHECK AND SAVE ============================================================
  # ==============================================================================
  
  # 3.1 Check --------------------------------------------------------------------
  
  # Check number of sites
  if(n_start != length(unique(ds.agg$site_code))){stop('After aggregating to site level number of sites does not match original number of sites')}

  # Check NAs
  # There will be some NAs introduced for: ccdc_nextSegmentStart, ccdc_nextSegmentEnd, ccdc_previousSegmentStart, ccdc_previousSegmentEnd
  # That is okay, this is where all plots within the site were 'no segment'
  # For full synthesis dataset ('all') NAs for plot size and biomass are okay, these are PFTs that were not recorded
  colSums(is.na(ds.agg))

  # 3.2 Save --------------------------------------------------------------------
  
  write.csv(ds.agg, paste0(out_dir, 'ds_', ds_name, '_site_level_', version, '.csv'), row.names = FALSE)

}
