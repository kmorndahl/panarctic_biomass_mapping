################################################################################
################################################################################

# DESCRIPTION:
# 
# This script creates the Pan Arctic region of interest that is used to 
# define the study bounds
#   - Cleans up existing Pan Arctic shapefile
#    - Repairs geometry
#    - Fills holes and slivers

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

# NOTES:

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

library(sf)
library(tidyverse)
library(smoothr)
library(units)
library(nngeo)

# 1.2 Parameters ---------------------------------------------------------------

coast_buffer_dist = 175000 # 175 km
tidy_threshold = units::set_units(1, km^2)
gap_threshold = units::set_units(1000, km^2)
icesheet_tidy_threshold = units::set_units(100, km^2)

in_dir = 'data/000_gis/'
out_dir = 'data/00_roi/'

# 1.3 Read in data -------------------------------------------------------------

greenland = sf::st_read(paste0(in_dir, 'GRL_adm0.shp'))
ocean = sf::st_read(paste0(in_dir, 'ne_10m_ocean.shp'))
arctic = sf::st_read(paste0(in_dir, 'arctic_oroarctic_laea.shp'))
greenland_icesheet = sf::st_read(paste0(in_dir, 'Greenland_Basins_PS_v1.4.2.shp'))

greenland = greenland %>% dplyr::select(geometry)
ocean = ocean %>% dplyr::select(geometry)
arctic = arctic %>% dplyr::select(geometry)
greenland_icesheet = greenland_icesheet %>% dplyr::select(geometry)

# ==============================================================================
# 2. TIDY GEOMETRIES ===========================================================
# ==============================================================================

# 2.1 Pan Arctic ---------------------------------------------------------------

arctic = sf::st_make_valid(arctic) # Repair geometry
arctic = smoothr::fill_holes(arctic, tidy_threshold) # Remove holes
arctic = sf::st_buffer(arctic, 1) # Buffer by 1 m to remove small slivers that are not technically holes

# 2.2 Greenland icesheet -------------------------------------------------------

greenland_icesheet = sf::st_make_valid(greenland_icesheet) # Repair geometry
greenland_icesheet = sf::st_union(greenland_icesheet) # Dissolve
greenland_icesheet = sf::st_buffer(greenland_icesheet, 50) # Buffer to remove small slivers lingering after dissolve
greenland_icesheet = sf::st_buffer(greenland_icesheet, -50) # Reverse buffer to return to original

# ==============================================================================
# 3. CREATE PAN ARCTIC ROI =====================================================
# ==============================================================================

# 3.1 Reproject ----------------------------------------------------------------

greenland = sf::st_transform(greenland, st_crs(arctic))
greenland_icesheet = sf::st_transform(greenland_icesheet, st_crs(arctic))
ocean = sf::st_transform(ocean, st_crs(arctic))

# 3.2 Fix Greenland -------------------------------------------------------------

# Union original arctic polygon + better defined Greenland boundary
arctic_grnlnd = sf::st_union(arctic, greenland)

# 3.3 Buffer -------------------------------------------------------------------

bufffer_shp_name = paste0('arctic_oroarctic_grnlnd_buf', coast_buffer_dist/1000, 'km_laea.shp')
arctic_grnlnd_buffer = sf::st_buffer(arctic_grnlnd, coast_buffer_dist)

# 3.4 Reduce buffer to coastline only ------------------------------------------

# Clip ocean layer to buffered arctic polygon to create coastal buffer
coast_buffer = sf::st_intersection(ocean, arctic_grnlnd_buffer)

# Join arctic polygon (NOT buffered) with coastal buffer
arctic_grnlnd_coast_buffer = sf::st_union(arctic_grnlnd, coast_buffer)

# Tidy
arctic_grnlnd_coast_buffer = arctic_grnlnd_coast_buffer %>% select(geometry)

# Fill holes 1
arctic_grnlnd_coast_buffer_fill_holes = smoothr::fill_holes(arctic_grnlnd_coast_buffer, gap_threshold)
arctic_grnlnd_coast_buffer_fill_holes = sf::st_make_valid(arctic_grnlnd_coast_buffer_fill_holes) # Repair geometry

# Fill holes 2
# (Need to perform twice in a row to thoroughly remove holes)
arctic_grnlnd_coast_buffer_fill_holes_2 = smoothr::fill_holes(arctic_grnlnd_coast_buffer_fill_holes, gap_threshold)
arctic_grnlnd_coast_buffer_fill_holes_2 = sf::st_make_valid(arctic_grnlnd_coast_buffer_fill_holes_2) # Repair geometry

# Fix slivers
arctic_grnlnd_coast_buffer_fix_slivers = sf::st_buffer(arctic_grnlnd_coast_buffer_fill_holes_2, 1) # Buffer to remove small slivers
arctic_grnlnd_coast_buffer_fix_slivers = sf::st_buffer(arctic_grnlnd_coast_buffer_fix_slivers, -1) # Reverse buffer to return to original

# 3.5 Remove Greenland icesheet -----------------------------------------------

# Remove Greenland icesheet
arctic_grnlnd_coast_buffer_rm_icesheet = sf::st_difference(arctic_grnlnd_coast_buffer, greenland_icesheet)

# Remove holes and drop crumbs
arctic_grnlnd_coast_buffer_rm_icesheet = smoothr::fill_holes(arctic_grnlnd_coast_buffer_rm_icesheet, icesheet_tidy_threshold)
arctic_grnlnd_coast_buffer_rm_icesheet = smoothr::drop_crumbs(arctic_grnlnd_coast_buffer_rm_icesheet, icesheet_tidy_threshold)

# ==============================================================================
# 4. SAVE FINAL SHAPEFILES =====================================================
# ==============================================================================

st_write(arctic_grnlnd, paste0(in_dir, 'arctic_oroarctic_grnlnd_laea.shp'))

st_write(arctic_grnlnd_coast_buffer_fix_slivers, paste0(out_dir, 'arctic_oroarctic_coast_buffer_laea.shp'))

st_write(arctic_grnlnd_coast_buffer_rm_icesheet, paste0(out_dir, 'arctic_oroarctic_polygon.shp'))
