################################################################################
################################################################################

# DESCRIPTION:

# This script creates figure of study area with bioclimate zones and field data

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

# 1.2 Parameters ---------------------------------------------------------------

pt_size = 5
cb_palette = c("#F5C710", "#21918c", "#942c80", "#117733") # Viridis colors

gis_dir = 'data/000_gis/'
data_dir = 'output/06_model_ready_data/'
out_dir = 'output/17_figures/'

# 1.3 Read in data -------------------------------------------------------------

deg45 = st_read(paste0(gis_dir, '45N.shp'))
zones = st_read(paste0(gis_dir, 'arctic_three_lat_zones_laea.shp'))
arctic = st_read(paste0(gis_dir, 'arctic_oroarctic_grnlnd_laea.shp'))
land = st_read(paste0(gis_dir, 'land_45n_laea.shp'))
land_all = st_read(paste0(gis_dir, 'ne_50m_land.shp'))
data = read.csv(paste0(data_dir, 'ds_woody_v20240508.csv')) %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# 1.4 Establish circular border ------------------------------------------------

north_pole = data.frame(lat = -760000, long = 0)
north_pole= st_as_sf(north_pole, coords = c("long", "lat"), crs = crs(land))
circle = st_buffer(north_pole, 4300000)
box = st_as_sfc(st_bbox(circle))
circular_frame = st_difference(box, circle)

# ==============================================================================
# 2. PLOT ======================================================================
# ==============================================================================

plt = ggplot() +
  geom_sf(data = land, fill = 'grey90')+
  geom_sf(data = zones, aes(fill = FIRST_zone))+
  geom_sf(data = circular_frame, color = 'white', fill = 'white')+
  geom_sf(data = data[data$dataset_id == 'nonveg',], fill = 'white', color = 'black', pch=21, size = pt_size)+
  geom_sf(data = data[data$dataset_id != 'nonveg',], fill = 'black', color = 'white', pch=21, size = pt_size)+
  coord_sf(ylim = c(st_bbox(circular_frame)$ymin + 400000, st_bbox(circular_frame)$ymax - 400000),
           xlim = c(st_bbox(circular_frame)$xmin + 400000, st_bbox(circular_frame)$xmax - 400000),
           default_crs = crs(circular_frame),
           lims_method = 'box')+
  scale_fill_manual(values = cb_palette)+
  scale_colour_manual(values = cb_palette)+
  theme_minimal(base_size = 65)+
  theme(legend.position = 'top', 
        legend.key.width = unit(1.7, "cm"), 
        legend.key.height = unit(1.7, "cm"),
        axis.text = element_text(size = 40),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(col = '', fill = '')
plt

# ==============================================================================
# 3. SAVE ======================================================================
# ==============================================================================

ggsave(
  paste0(out_dir, 'zones_map.jpg'),
  plt,
  width = 40,
  height = 40,
  units = 'cm',
  bg = 'white',
  dpi = 600
)