################################################################################
################################################################################

# DESCRIPTION:

# This script creates 3D hillshade of biomass map

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

library(tidyverse)
library(raster)
library(terra)
library(cmocean)
library(rayshader)
library(magick)
library(ggtern)
library(sf)

# 1.2 Parameters ---------------------------------------------------------------

application = 'figure' # Choose 'figure' or 'graphical_abstract'

if(application == 'graphical_abstract'){
  zscale = 4 # Ratio of DEM vertical resolution (1 m) and horizontal resolution (2 m - 90 m depending on DEM source)
  theta = 45
  phi = 45
  fov = 0
  transect_lwd = 10
  roi = 'ns'
}else if(application == 'figure'){
  zscale = 2 # Ratio of DEM vertical resolution (1 m) and horizontal resolution (2 m - 90 m depending on DEM source)
  theta = 330
  phi = 45
  fov = 0
  transect_lwd = 5
  roi = 'ural'
}else{
  stop('Application not recognized, please choose "figure" or "graphical_abstract"')
}

in_dir = 'output/15_gee_output/topography/'
out_dir = 'output/17_figures/'

# 1.3 Read in data -------------------------------------------------------------

biomass = terra::rast(paste0(in_dir, 'hillshade_', roi, '_biomass.tif'))
dem = terra::rast(paste0(in_dir, 'hillshade_', roi, '_arcticdem.tif'))
transect = read_sf(paste0(in_dir, 'biomass_elev_transect_fc_', roi, '.shp'))
water = read_sf(paste0(in_dir, 'NHDArea.shp'))

# ==============================================================================
# 2. FORMAT DATA ===============================================================
# ==============================================================================

# 2.1 Define extent ------------------------------------------------------------

if(application == 'graphical_abstract'){
  scale_xmin = 0
  scale_xmax = 0
  scale_ymin = 0
  scale_ymax = 0
}else if(application == 'figure'){
  scale_xmin = 7000
  scale_xmax = 8500
  scale_ymin = 7000
  scale_ymax = 8500
}else{
  stop('Application not recognized, please choose "figure" or "graphical_abstract"')
}

# Assign boundary coordinates
lon = c(st_bbox(biomass)$xmin[[1]]+scale_xmin, st_bbox(biomass)$xmax[[1]]-scale_xmax)
lat = c(st_bbox(biomass)$ymin[[1]]+scale_ymin, st_bbox(biomass)$ymax[[1]]-scale_ymax)
coord_df = data.frame(lon, lat)

# Convert to spatial object
bbox = coord_df %>% 
  st_as_sf(coords = c("lon", "lat"), 
           crs = crs(biomass)) %>% 
  st_bbox() %>%
  st_as_sfc()

# 2.2 Convert transect from points to line -------------------------------------

transect$id = 1
transect = transect %>%
  group_by(id) %>%
  dplyr::summarize() %>%
  st_cast("LINESTRING") 

# 2.3 Reproject and clip -------------------------------------------------------

biomass = crop(biomass, bbox)
dem = crop(dem, bbox)
transect = st_transform(transect, crs(biomass))
transect = st_intersection(transect, bbox)

if(application == 'graphical_abstract'){
  water = st_transform(water, crs(biomass))
  water = st_intersection(water, bbox)
}

# 2.4 Convert rasters to RGB ---------------------------------------------------
# Necessary in order to use them as overlays

# Cap biomass at 1,500 g/m2 to match scale in figures
if(application == 'figure'){
  biomass[biomass$plant_biomass_gm2>1500]  = 1500
}

# Convert rasters to matrices
biomass_df = as.data.frame(biomass, xy = TRUE) 
biomass_mat = raster_to_matrix(biomass)
dem_df = as.data.frame(dem, xy = TRUE) 
dem_mat = raster_to_matrix(dem)

# Assign biomass values to palette values
biomass_min = min(biomass_df$plant_biomass_gm2)
biomass_max = max(biomass_df$plant_biomass_gm2)
coltb = data.frame(value = 0:biomass_max, col = cmocean('speed')(biomass_max + 1))

# Assign to raster color table
terra::coltab(biomass) = coltb 

# Extract new color table as dataframe
biomass_coltb = data.frame(terra::coltab(biomass)) 

# Convert to RGB hex codes
hex_code = ggtern::rgb2hex(
  r = biomass_coltb[,2],
  g = biomass_coltb[,3],
  b = biomass_coltb[,4]
)

# Assign RGB colors to raster
cols = hex_code
from = 0:biomass_max
to = t(col2rgb(cols))
biomass = na.omit(biomass)
biomass_rgb = terra::subst(
  biomass,
  from = from,
  to = to,
  names = cols
)

# Save and load RGB raster
# Necessary to save and read in again for it to work properly as an overlay
out_path = paste0(out_dir, 'hillshade_', roi, '_biomass_rgb.png')
terra::writeRaster(
  biomass_rgb,
  out_path,
  overwrite = T,
  NAflag = 255
)
img = png::readPNG(out_path)

# 2.5 Create overlays ----------------------------------------------------------

# Create water overlay
if(application == 'graphical_abstract'){
  overlay = generate_polygon_overlay(
    geometry = water,
    extent = st_bbox(biomass),
    heightmap = dem_mat,
    linecolor = NA,
    palette = "lightblue",
    linewidth = 0
  )
}

# Create transect overlay
if(application == 'figure'){
  overlay = generate_line_overlay(
    geometry = transect,
    extent = st_bbox(biomass),
    heightmap = dem_mat,
    color = "magenta",
    linewidth = transect_lwd,
    lty = 
  )
}

# ==============================================================================
# 3. PLOT ======================================================================
# ==============================================================================

dem_mat %>%
  sphere_shade() %>%
  add_shadow(ray_shade(dem_mat, zscale = zscale), 0.5) %>%
  add_overlay(overlay = img, alphalayer = 0.8) %>%
  add_overlay(overlay) %>%
  plot_3d(dem_mat, zscale = zscale, theta = theta, phi = phi, fov = fov)

Sys.sleep(0.2)
render_snapshot(paste0(out_dir, 'hillshade_3d_', roi))
