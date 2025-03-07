################################################################################
################################################################################

# DESCRIPTION:

# This script produces a plot of the spatial cross-validation blocking

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

# NOTES:

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

library(spatialsample)
library(sf)

# 1.2 Parameters ---------------------------------------------------------------

set.seed(1908)

version = 'v20240508'
biomass_dir = 'output/06_model_ready_data/'
gis_dir = 'data/000_gis/'
out_fig_dir = 'output/07_exploratory_analysis/spatial_blocking/'

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

biomass_df = read.csv(paste0(biomass_dir, 'ds_woody_', version, '.csv'))
arctic = sf::st_read(paste0(gis_dir, 'arctic_oroarctic_clean_laea.shp'))
land = sf::st_read(paste0(gis_dir, 'land_45n_laea.shp'))

# ==============================================================================
# 2. TIDY ======================================================================
# ==============================================================================

# 2.1 Finalize biomass data ----------------------------------------------------

# Convert data to sf
biomass_sf = sf::st_as_sf(biomass_df, coords = c("longitude", "latitude"), crs = 4326)

# Transform biomass data
biomass_sf = st_transform(biomass_sf, 3571)

# 2.2 Establish circular border ------------------------------------------------

north_pole = data.frame(lat = -760000, long = 0)
north_pole= st_as_sf(north_pole, coords = c("long", "lat"), crs = crs(land))
circle = st_buffer(north_pole, 4300000)
box = st_as_sfc(st_bbox(circle))
box = st_buffer(box, 900000)
circular_frame = st_difference(box, circle)

# ==============================================================================
# 3. CREATE FOLDS ==============================================================
# ==============================================================================

# Create spatial cross-validation folds
block_folds = spatialsample::spatial_block_cv(biomass_sf, v = NULL, cellsize = 100*1000)

# ==============================================================================
# 4. PLOT ======================================================================
# ==============================================================================

b = spatialsample::autoplot(block_folds, alpha = 0.75, size = 1, show_grid = TRUE) + 
  geom_sf(data = circular_frame, color = 'white', fill = 'white')+
  geom_sf(data = arctic, fill = 'grey80', lwd = 0.1)+
  geom_sf(data = land, fill = 'grey90', lwd = 0.1)+
  coord_sf(ylim = c(st_bbox(circle)$ymin + 400000, st_bbox(circle)$ymax - 400000),
           xlim = c(st_bbox(circle)$xmin + 400000, st_bbox(circle)$xmax - 400000),
           default_crs = st_crs(circle),
           lims_method = 'box')+
  theme_minimal(base_size = 30)+
  theme(legend.position = 'none', 
        axis.text = element_text(size = 30),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(col = '', fill = '')

b$layers = rev(b$layers)
b

ggsave(
  paste0(out_fig_dir, '/blocking_100km_gray.png'),
  b,
  bg = 'white',
  dpi = 300
)

