################################################################################
################################################################################

# DESCRIPTION:

# This script creates final dataset with plot locations and dates for use in GEE

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
require(tidyverse)
require(data.table)

# 1.2 Parameters ---------------------------------------------------------------

version = 'v20240215'
synth_dir = 'data/01_synthesis_dataset/'
dir = 'data/02_gee_ready/'

# 1.3 Read in data -------------------------------------------------------------

ds.synth.all = fread(paste0(synth_dir, 'the_arctic_plant_aboveground_biomass_synthesis_dataset_v1.2.csv'))
ds.plots = fread(paste0(dir, 'arctic_tundra_biomass_synthesis_plots_', version, '.csv'))

# SUBSET DATA AND REMOVE DUPLICATE PLOTS ======================================================

ds.synth.all.dates = ds.synth.all %>% 
  select('site_code', 'plot_code', 'year', 'month', 'day') %>%
  distinct(plot_code, .keep_all = TRUE) %>%
  as.data.table()

# ==============================================================================
# 2. TIDY AND COMBINE ==========================================================
# ==============================================================================

# 2.1 Add midsummer date for non-vegetated data --------------------------------

ds.nonveg = ds.plots[ds.plots$dataset_id == 'nonveg', list(site_code, plot_code, year)]
ds.nonveg$month = 7 # Choose midsummer
ds.nonveg$day = 31 # Choose midsummer

# 2.2 Join user created non-vegetated data and harvest plot data ---------------

ds.synth.all.dates = rbind(ds.synth.all.dates, ds.nonveg)

# ==============================================================================
# 3. SAVE ======================================================================
# ==============================================================================

write.csv(ds.synth.all.dates, paste0(dir, 'arctic_tundra_biomass_synthesis_dataset_plot_dates_', version, '.csv'), row.names = FALSE)