################################################################################
################################################################################

# DESCRIPTION:

# This script adds combines aboveground vegetation biomass synthesis dataset
# data with user created non-vegetated observations and tidies for use in GEE

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
require(sf)
require(tidyverse)
require(data.table)

version = 'v20240215'

# 1.2 Parameters ---------------------------------------------------------------

synth_dir = 'data/01_synthesis_dataset/'
out_dir = 'data/02_gee_ready/'

# 1.3 Read in data -------------------------------------------------------------

ds.synth.plots = fread(paste0(synth_dir, 'arctic_tundra_biomass_synthesis_dataset_plot_locations.csv'))
ds.nonveg = st_read(paste0(synth_dir, 'barren_built_snow_ice_water_random_pts_scale30m.shp'))

# ==============================================================================
# 2. TIDY AND COMBINE ==========================================================
# ==============================================================================

# 2.1 Tidy non-vegetated data --------------------------------------------------

# Label non-vegetated types
ds.nonveg$Map = as.character(ds.nonveg$Map)
ds.nonveg = ds.nonveg %>% dplyr::mutate(lc = dplyr::recode(Map, '50'='built', '60'='barrens', '70'='snowice', '80'='water')) %>% select(-Map)

# Select valid non-vegetated points
ds.nonveg = ds.nonveg[ds.nonveg$acceptable == 'yes',]

# Write lat/long to fields
ds.nonveg = cbind(ds.nonveg %>% st_drop_geometry(), # Drop geometry
                  data.frame(st_coordinates(ds.nonveg))) %>% # Get geometry as coordinates
                  select(-acceptable) # Remove 'quality' column

# Populate fields to match plot data
ds.nonveg.tidy = data.frame(dataset_id = rep('nonveg', nrow(ds.nonveg)),
                            citation_short = rep('nonveg', nrow(ds.nonveg)),
                            site_code = paste0('nonveg_', ds.nonveg$id),
                            plot_code = paste0(ds.nonveg$lc, '_', ds.nonveg$id),
                            coord_type = rep('plot', nrow(ds.nonveg)),
                            year = 2021, # ESA WorldCover data is from 2021, QGIS Bing and Google satellite imagery is presumably fairly recent
                            latitude = ds.nonveg$Y,
                            longitude = ds.nonveg$X)

# 2.2 Join user created non-vegetated data and harvest plot data ---------------

ds.plots = rbind(ds.synth.plots, ds.nonveg.tidy)

# ==============================================================================
# 3. SAVE ======================================================================
# ==============================================================================

write.csv(ds.plots, paste0(out_dir, 'arctic_tundra_biomass_synthesis_plots_', version, '.csv'), row.names = FALSE)
