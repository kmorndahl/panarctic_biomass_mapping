################################################################################
################################################################################

# DESCRIPTION:

# This script creates figure of study area with bioclimate zones and map tile grid

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

# NOTES:

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)
library(raster)
library(fasterize)
library(smoothr)
library(units)

# 1.2 Parameters -----------------------------------------------------------------

cb_palette = c("#F5C710", "#21918c", "#942c80", "#117733")

in_dir = 'data/000_gis/'
out_dir = 'output/18_final_maps/'

# 1.3 Read in data -------------------------------------------------------------

deg45 = st_read(paste0(in_dir, '45N.shp'))
arctic = st_read(paste0(in_dir, 'arctic_oroarctic_grnlnd_laea.shp'))
zones = st_read(paste0(in_dir, 'arctic_three_lat_zones_laea.shp'))
land = st_read(paste0(in_dir, 'land_45n_laea.shp'))
land_all = st_read(paste0(in_dir, 'ne_50m_land.shp'))
grid = st_read(paste0(out_dir, 'tile_index.shp'))


################################################################################
################################################################################

# ==============================================================================
# 2. PLOT ======================================================================
# ==============================================================================

# 2.1 Establish circular border ------------------------------------------------

north_pole = data.frame(lat = -760000, long = 0)
north_pole= st_as_sf(north_pole, coords = c("long", "lat"), crs = crs(land))
circle = st_buffer(north_pole, 4300000)
box = st_as_sfc(st_bbox(circle))
circular_frame = st_difference(box, circle)

# 2.2 Plot -------------------------------------------------------------

plt = ggplot() +
  geom_sf(data = land, fill = 'grey90')+
  geom_sf(data = zones, aes(fill = FIRST_zone))+
  geom_sf(data = grid, fill = NA, linewidth = 1.5)+
  geom_sf_text(data = grid, aes(label = tile_id), size = 16) +
  scale_fill_manual(values = cb_palette)+
  scale_colour_manual(values = cb_palette)+
  theme_minimal(base_size = 50)+
  theme(legend.position = 'top', 
        legend.key.width = unit(1.7, "cm"), 
        legend.key.height = unit(1.7, "cm"),
        axis.text = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y  =element_blank(),
        panel.grid = element_line(colour = "grey90", linewidth = 0.5))+
  labs(col = '', fill = '', x = '', y = '')
plt


################################################################################
################################################################################

# ==============================================================================
# 3. SAVE ======================================================================
# ==============================================================================

ggsave(
  paste0(out_dir, 'tile_index_map.jpg'),
  plt,
  width = 40,
  height = 40,
  units = 'cm',
  bg = 'white',
  dpi = 600
)