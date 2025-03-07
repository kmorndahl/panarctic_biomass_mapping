################################################################################
################################################################################

# DESCRIPTION:
# 
# This script creates the final bioclimate zone shapefiles
#   - Adds better constrained Greenland boundary
#   - Adds subarctic zone

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
library(smoothr)
library(units)

# 1.2 Parameters -----------------------------------------------------------------

in_dir = 'data/000_gis/'
out_dir = 'data/00_rois/'

# 1.3 Read in data -------------------------------------------------------------

zones = sf::st_read(paste0(in_dir, 'arctic_three_lat_zones_laea.shp'))
arctic = sf::st_read(paste0(in_dir, 'arctic_oroarctic_grnlnd_laea.shp'))
land = sf::st_read(paste0(in_dir, 'land_45n_laea.shp'))
greenland = sf::st_read(paste0(in_dir, 'GRL_adm0.shp') %>% dplyr::select(geometry)

# ==============================================================================
# 2. TIDY ======================================================================
# ==============================================================================

# 2.1 Fix sliver ---------------------------------------------------------------

zones = sf::st_buffer(zones, 1) # Buffer to remove small slivers lingering after dissolve
zones = sf::st_buffer(zones, -1) # Reverse buffer to return to original

# 2.2 Fix Greenland ------------------------------------------------------------

# Transform
greenland = sf::st_transform(greenland, st_crs(zones))

# Get High Arctic
high_arctic = zones[1,]

# Union High Arctic + better defined Greenland boundary
high_arctic_grnlnd = sf::st_union(high_arctic, greenland)

# Get Low Arctic
low_arctic = zones[2,]

# Difference to remove areas of the new High Arctic + Greenland that are
# actually part of the Low Arctic
high_arctic_grnlnd_final = sf::st_difference(high_arctic_grnlnd, low_arctic)
high_arctic_grnlnd_final = high_arctic_grnlnd_final %>% dplyr::select(c(dsl, FIRST_zone, geometry))

# Replace High Arctic
zones[1,] = high_arctic_grnlnd_final

# Save 
st_write(zones, paste0(out_dir, 'zones_grnlnd_laea.shp'))

# 2.3 Add subarctic -------------------------------------------------------------

# Clip land by arctic, the remainder is the "subarctic"
land_clipped = st_difference(land, arctic)

# Clean up "crumbs" due to mismatched coastline
subarctic = smoothr::drop_crumbs(land_clipped, units::set_units(11000000000, m^2))

# Add subarctic to arctic zones
subarctic = subarctic %>% rename(dsl = dslv, FIRST_zone = FID)
subarctic$dsl = 4
subarctic$FIRST_zone = 'Sub Arctic'
zones = bind_rows(zones, subarctic)
st_write(zones, paste0(out_dir, 'zones_grnlnd_laea_subarctic.shp'))



