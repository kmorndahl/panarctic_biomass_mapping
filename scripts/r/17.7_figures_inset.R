################################################################################
################################################################################

# DESCRIPTION:

# This script creates inset figures showing case study locations

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

# 1.2 Parameters ---------------------------------------------------------------

fig = 'disturbance' # Choose 'disturbance', 'seasonal_refl', or 'topo_corr'
pt_size = 18

in_dir = 'data/000_gis/'
out_dir = 'output/17_figures/'

# 1.3 Read in data -------------------------------------------------------------

land = st_read(paste0(in_dir, 'land_45n_laea.shp'))

if(fig == 'disturbance'){
  fires = st_read(paste0(in_dir, 'jones_fires.shp'))
  thaw_slumps = st_read(paste0(in_dir, 'thaw_slumps.shp'))
  meade = st_centroid(fires)[2,] %>% select(geometry)
  ketik = st_centroid(fires)[1,] %>% select(geometry)
  horton = st_centroid(thaw_slumps[thaw_slumps$id == 716,]) %>% select(geometry)
  roi = bind_rows(meade, horton)
}else if(fig == 'seasonal_refl'){
  roi = data.frame(lat = c(69.9964), long = c(27.9818))
  roi = roi %>% st_as_sf(coords = c("long", "lat"), crs = 4326)
}else if(fig == 'topo_corr'){
  roi = data.frame(lat = c(68.2249), long = c(138.1972))
  roi = roi %>% st_as_sf(coords = c("long", "lat"), crs = 4326)
}else{
  stop('Figure type not recognized')
}

# 1.4 Establish circular border ------------------------------------------------

north_pole = data.frame(lat = -760000, long = 0)
north_pole= st_as_sf(north_pole, coords = c("long", "lat"), crs = st_crs(land))
circle = st_buffer(north_pole, 4300000)
box = st_as_sfc(st_bbox(circle))
circular_frame = st_difference(box, circle)

# ==============================================================================
# 2. PLOT ======================================================================
# ==============================================================================

plt = ggplot() +
  geom_sf(data = land, fill = 'grey90')+
  geom_sf(data = roi, fill = 'magenta', color = 'black', pch=21, size = 18)+
  theme_minimal(base_size = 60)+
  theme(legend.position = NULL, 
        axis.text = element_text(size = 40),
        line = element_blank(), 
        rect = element_blank(), 
        text = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank())+
  labs(col = '', fill = '', x = '', y = '')

# ==============================================================================
# 3. SAVE ======================================================================
# ==============================================================================

ggsave(
  paste0(out_dir, fig, '_inset_map.jpg'),
  plt,
  width = 40,
  height = 40,
  units = 'cm',
  bg = 'white',
  dpi = 600
)