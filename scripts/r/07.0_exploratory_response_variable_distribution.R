################################################################################
################################################################################

# DESCRIPTION:

# This script produces histograms of response variable (biomass)

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
require(ggridges)

# 1.2 Parameters ---------------------------------------------------------------

in_dir = 'output/06_model_ready_data/'
version = 'v20240508'
data_subset = 'presence_only' # Choose 'all_data' or 'presence_only'
out_fig_dir = 'output/07_exploratory_analysis/histograms/'

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

if(data_subset == 'presence_only'){
  ds.woody = ds.woody[ds.woody$biomass_density_gm2 > 0,]
  ds.total = ds.total[ds.total$biomass_density_gm2 > 0,]
}

# 2.3 Histograms ------------------------------------------------

cb_palette = c("#88CCEE", "#CC6677", "#888888", "#117733")

hist_total = ggplot(ds.total, aes(x =  biomass_density_gm2))+
  geom_histogram()+
  labs(y = 'Count', x = bquote('Total Plant Biomass'~(g~m^-2)))+
  theme_minimal(base_size = 30)
hist_total

hist_woody = ggplot(ds.woody, aes(x =  biomass_density_gm2))+
  geom_histogram()+
  labs(y = 'Count', x = bquote('Woody Plant Biomass'~(g~m^-2)))+
  theme_minimal(base_size = 30)
hist_woody

dens_total = ggplot(ds.total, aes(x = biomass_density_gm2, y = zone_3, fill = zone_3))+
  geom_density_ridges()+
  scale_fill_manual(values = cb_palette)+
  labs(y = 'Zone', x = bquote('Total Plant Biomass'~(g~m^-2)))+
  theme_minimal(base_size = 30)+
  theme(legend.position = "none")
dens_total

dens_woody = ggplot(ds.woody, aes(x = biomass_density_gm2, y = zone_3, fill = zone_3))+
  geom_density_ridges()+
  scale_fill_manual(values = cb_palette)+
  labs(y = 'Zone', x = bquote('Woody Plant Biomass'~(g~m^-2)))+
  theme_minimal(base_size = 30)+
  theme(legend.position = "none")
dens_woody

# ==============================================================================
# 3. SAVE ======================================================================
# ==============================================================================

ggsave(
  paste0(out_fig_dir, 'histogram_total_',  data_subset, '.png'),
  hist_total,
  width = 40,
  height = 30,
  units = 'cm',
  bg = 'white',
  dpi = 600
)

ggsave(
  paste0(out_fig_dir, 'histogram_woody_',  data_subset, '.png'),
  hist_woody,
  width = 40,
  height = 30,
  units = 'cm',
  bg = 'white',
  dpi = 600
)

ggsave(
  paste0(out_fig_dir, 'density_zone_total_',  data_subset, '.png'),
  dens_total,
  width = 40,
  height = 30,
  units = 'cm',
  bg = 'white',
  dpi = 600
)

ggsave(
  paste0(out_fig_dir, 'density_zone_woody_',  data_subset, '.png'),
  dens_woody,
  width = 40,
  height = 30,
  units = 'cm',
  bg = 'white',
  dpi = 600
)