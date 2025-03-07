################################################################################
################################################################################

# DESCRIPTION:

# This script:
#   - calculates annual CCDC metrics (e.g. annual median, change between seasons)
#   - calculates tree cover using Hansen data
#   - calculates CCDC segment metadata
#   - join biomass harvest data with predictor data

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
require(R.utils)
require(lubridate)

# 1.2 Parameters ---------------------------------------------------------------

rep_version = 'v20240215'
predictor_version = 'topocorr_v20240207'
out_version = 'v20240319'
ccdc_start_year = 1984
ccdc_end_year = 2023

rep_dir = 'output/03_representativeness_analysis/joined/'
predictors_dir = 'output/04_predictor_extraction/'
out_dir = 'output/04_predictor_extraction/'

ds_name_list = c('total', 'woody')

# 1.3 Functions ----------------------------------------------------------------

find_mode = function(x) {
  u = unique(na.omit(x))
  tab = tabulate(match(x, u))
  u[tab == max(tab)][1] # If there are two modes, grab the first one
}

std.error = function(x) sd(na.omit(x))/sqrt(length(na.omit(x)))

# 1.4 Read in data -------------------------------------------------------------

for(ds_name in ds_name_list){
  
  # Read in synthesis database with biomass harvest information and representativeness data
  ds.synth = fread(paste0(rep_dir, 'arctic_tundra_', ds_name, '_biomass_synthesis_representativeness_joined_', rep_version, '.csv'))

  # Read in datasets from GEE with predictor information
  ds.predictors.plots = fread(paste0(predictors_dir, 'arctic_tundra_biomass_synthesis_plots_predictors_', ccdc_start_year, '_', ccdc_end_year, '_closest_', predictor_version, '.csv'))
  ds.predictors.sites = fread(paste0(predictors_dir, 'arctic_tundra_biomass_synthesis_sites_predictors_', ccdc_start_year, '_', ccdc_end_year, '_closest_', predictor_version, '.csv'))
  
  # ==============================================================================
  # 2. PREDICTOR DATA ============================================================
  # ==============================================================================
  
  # 2.1 Tidy -----------------------------------------------------------------
  
  # For plot level data, assign distance to pixel centroid as 0 (not exactly correct but just a placeholder)
  ds.predictors.plots$distance = 0 
  
  # For site level data, change plot code to match site code to indicate these predictor data are only available for one location per site
  ds.predictors.sites$plot_code = ds.predictors.sites$site_code
  
  # Remove problematic columns
  ds.predictors.plots = ds.predictors.plots %>% select(-c('.geo', 'system:index')) %>% as.data.table()
  ds.predictors.sites = ds.predictors.sites %>% select(-c('.geo', 'system:index')) %>% as.data.table()
  
  # Combine predictor data
  ds.predictors = rbind(ds.predictors.plots, ds.predictors.sites)
  
  # 2.2 Calculate CCDC annual metrics --------------------------------------------
  
  spectral = names(ds.predictors.plots)[grepl('spectral_', names(ds.predictors.plots))]
  spectral_bands = gsub('spectral_', '', spectral)
  spectral_bands = unique(gsub('_startSnowfree|_earlySummer|_peakSummer|_lateSummer|_endSnowfree', '', spectral_bands))
  
  for(band in spectral_bands){
    
    print(paste0('The band is: ', band))
    
    # Get spectral band names
    start_snowfree = paste0('spectral_', band, '_startSnowfree')
    early_summer = paste0('spectral_', band, '_earlySummer')
    peak_summer = paste0('spectral_', band, '_peakSummer')
    late_summer = paste0('spectral_', band, '_lateSummer')
    end_snowfree = paste0('spectral_', band, '_endSnowfree')
    
    # Get seasonal DOY names
    start_snowfree_doy = 'doy_startSnowfree'
    early_summer_doy = 'doy_earlySummer'
    peak_summer_doy = 'doy_peakSummer'
    late_summer_doy = 'doy_lateSummer'
    end_snowfree_doy = 'doy_endSnowfree'
    
    # Assign output band names
    annualMean = paste0('spectral_', band, '_annualMean')
    annualMedian = paste0('spectral_', band, '_annualMedian')
    annualMin = paste0('spectral_', band, '_annualMin')
    annualMax = paste0('spectral_', band, '_annualMax')
    annualRange = paste0('spectral_', band, '_annualRange')
    changeSSES = paste0('spectral_', band, '_changeSSES')
    changeESPS = paste0('spectral_', band, '_changeESPS')
    changePSLS = paste0('spectral_', band, '_changePSLS')
    changeLSES = paste0('spectral_', band, '_changeLSES')
    
    # Calculate annual summaries and seasonal rates of change
    # NOTE: take absolute value of day change for rates
    #   - for some barrens, DOYs are inconsistent i.e. sometimes peak summer is actually a day or so later than late summer
    #   - not a long term problem, usually only a day difference, likely due to prevalence of snow/ice, lack of data
    
    ds.predictors = ds.predictors %>%
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
    
  }
  
  # Replace NaN and Inf with zeros, this is where slopes have 0 denominators
  ds.predictors = ds.predictors %>% mutate(across(everything(), ~replace(.x, is.nan(.x), 0)))
  ds.predictors = ds.predictors %>% mutate(across(everything(), ~replace(.x, is.infinite(.x), 0)))
  
  # 2.3 Calculate tree loss ------------------------------------------------------
  
  ds.predictors$trees_cover = ds.predictors$trees_cover2000 
  
  # If the harvest year is after the loss year, zero cover is assigned
  ds.predictors$trees_cover[ds.predictors$year > (ds.predictors$trees_lossyear + 2000)] = 0
  
  # Tidy tree predictor
  # NAs indicate data at very high latitudes, outside of range of Hansen tree cover map
  # For these observations we can safely assign zero tree cover
  ds.predictors$trees_cover[is.na(ds.predictors$trees_cover)] = 0
  
  # Create binary tree cover column
  ds.predictors$trees_presence = ifelse(ds.predictors$trees_cover==0, 0, 1)
  
  # Account for gain up to 2012
  ds.predictors$trees_presence[(ds.predictors$trees_gain == 1) & ((ds.predictors$trees_lossyear + 2000) >= 2012)] = 1
  
  # Remove unnecessary columns
  ds.predictors = ds.predictors %>% select(!c(trees_cover2000, trees_loss, trees_gain, trees_lossyear))
  
  # 2.4 Assign CCDC segment metadata ---------------------------------------------
  
  # Get original number of plots and site pixels
  n_original_plots = nrow(ds.predictors)
  
  # Indicate which plots/sites are missing either previous or next segments
  ds.predictors$ccdc_previousSegmentStart[is.na(ds.predictors$ccdc_previousSegmentStart)] = 'no segment'
  ds.predictors$ccdc_previousSegmentEnd[is.na(ds.predictors$ccdc_previousSegmentEnd)] = 'no segment'
  ds.predictors$ccdc_nextSegmentStart[is.na(ds.predictors$ccdc_nextSegmentStart)] = 'no segment'
  ds.predictors$ccdc_nextSegmentEnd[is.na(ds.predictors$ccdc_nextSegmentEnd)] = 'no segment'
  
  # Indicate which plots/sites are using a segment that exactly matches the harvest date
  # Matt Macander: "if the sample date occurs within a segment the previous and next segment will be the same; but if the sample date occurs between segments the previous and next segment will be different"
  ds.predictors$ccdc_segmentType[ds.predictors$ccdc_previousSegmentStart == ds.predictors$ccdc_nextSegmentStart] = 'current'
  
  # If a plot/site is NOT using a segment that exactly matches the harvest date, indicate whether the previous or next segment is used
  ds.predictors$ccdc_segmentType[(ds.predictors$ccdc_previousSegmentStart != ds.predictors$ccdc_nextSegmentStart) & (ds.predictors$ccdc_tStart == ds.predictors$ccdc_previousSegmentStart)] = 'previous'
  ds.predictors$ccdc_segmentType[(ds.predictors$ccdc_previousSegmentStart != ds.predictors$ccdc_nextSegmentStart) & (ds.predictors$ccdc_tStart == ds.predictors$ccdc_nextSegmentStart)] = 'next'
  
  # ==============================================================================
  # 3. SYNTHESIS DATA ============================================================
  # ==============================================================================
  
  # Split synthesis data into plot and site level
  
  ds.synth.plots = ds.synth[ds.synth$coord_type == 'plot',]
  ds.synth.sites = ds.synth[ds.synth$coord_type == 'site',]

  # ==============================================================================
  # 4. COMBINE DATA ==============================================================
  # ==============================================================================
  
  # Identify column names duplicated across synthesis dataset and predictors dataset
  matching_cols = intersect(names(ds.synth), names(ds.predictors))
  
  # 4.1 Plot level ---------------------------------------------------------------
  
  ds.plots = merge(ds.synth.plots, ds.predictors %>% select(-matching_cols[matching_cols != "plot_code"]), by = c("plot_code"), all.x = TRUE) 
  
  # 4.2 Site level ---------------------------------------------------------------
  
  ds.sites = merge(ds.synth.sites, ds.predictors %>% select(-matching_cols[matching_cols != "site_code"]), by = c("site_code"), all.x = TRUE, allow.cartesian=TRUE) 

  # ==============================================================================
  # 5. TIDY DATA =================================================================
  # ==============================================================================
  
  # 5.1 Recombine plot and site level data ---------------------------------------
  
  ds = rbind(ds.plots, ds.sites)

  # 5.2 Derive additional phenological predictors --------------------------------
  
  # Add growing season length
  ds$doy_lengthSnowfree = ds$doy_endSnowfree - ds$doy_startSnowfree

  # Add summer length
  ds$doy_lengthSummer = ds$doy_lateSummer - ds$doy_earlySummer

  # Add day of year
  ds$day[is.na(ds$day)] = 15 # Where day is missing, assign middle of the month

  # Calculate harvest day-of-year
  ds$harvest_doy = yday(make_date(ds$year, ds$month, ds$day))

  # 5.3 Filter out plots with no spectral data -----------------------------------
  
  # Plots/sites with no CCDC fits and no spectral data: 
  #   - These are typically all nonveg sites in areas with a lot of snow/ice/clouds
  # Plots/sites with CCDC fits, but no spectral data:
  #   - This occurs when there is no segment for the harvest date + the closest segment is greater than extrapolateMaxDays=120 from the harvest date
  ds = ds[!is.na(ds$spectral_NDVI_peakSummer),]

  # 5.4 Manually fill world terrestrial ecosystem NAs ----------------------------
  
  # 14 = Polar Dry Sparsely or Non-vegetated on Plains
  ds$world_terrestrial_ecosystems[ds$site_code == 'USA.NorthSlope.HoweIsland-C12' & is.na(ds$world_terrestrial_ecosystems)] = 14
  ds$world_terrestrial_ecosystems[ds$site_code == 'USA.NorthSlope.HoweIsland-C13' & is.na(ds$world_terrestrial_ecosystems)] = 14

  # 5.5 Check NAs ----------------------------------------------------------------
  
  # For the full synthesis dataset ('all') NAs for biomass and plot area are okay -- these are unmeasured PFTs
  colnames(ds[ds$dataset_id != 'nonveg',])[apply(ds[ds$dataset_id != 'nonveg',], 2, anyNA)] 

  # ==============================================================================
  # 6. SAVE ======================================================================
  # ==============================================================================
  
  fwrite(ds, paste0(out_dir, 'arctic_tundra_', ds_name, '_biomass_synthesis_representativeness_predictors_', out_version, '.csv'))

}