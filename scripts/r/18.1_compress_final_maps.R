################################################################################
################################################################################

# DESCRIPTION:

# This script compresses final map rasters

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

library(terra)

# 1.2 Parameters ---------------------------------------------------------------

in_dir = '/agb/plant'

# 1.3 Read in data -------------------------------------------------------------

grid_lookup = read.csv('output/18_final_maps/0_grid_renaming.csv')
in_paths_agb = list.files(paste0('output/18_final_maps/agb/'), full.names = TRUE, recursive = TRUE)
in_paths_mask = list.files(paste0('output/18_final_maps/mask/'), full.names = TRUE, recursive = TRUE)
in_paths = c(in_paths_agb, in_paths_mask)

################################################################################
################################################################################

# ==============================================================================
# 2. COMPRESS AND SAVE =========================================================
# ==============================================================================

for(path in in_paths){
  
  print(paste0("Current raster: ", path))
  
  # Read in raster
  cog = rast(path)
  
  # Tidy path name
  out_path = gsub("plant_biomass", "plant_agb", path)
  out_path = gsub("p2_5", "p025", out_path)
  out_path = gsub("p50", "p500", out_path)
  out_path = gsub("p97_5", "p975", out_path)
  out_path = gsub("woody_plant", "woodyplant", out_path)
  
  # Get old grid cell name
  grid_cell = strsplit(strsplit(tail(strsplit(out_path, "_")[[1]], 1), "[.]")[[1]][1], "-")[[1]][2:3]
  grid_cell = paste(grid_cell, collapse = "-")
  
  # Get new grid cell name
  new_grid_cell = grid_lookup$new_name[grid_lookup$old_name == grid_cell]
  
  # Final tidying
  out_path = gsub(paste0("-", grid_cell), paste0("_", new_grid_cell), out_path)

  # Compress and write out
  writeRaster(cog, out_path, gdal=c("COMPRESS=LZW"), overwrite = FALSE)
  
  # Delete old raster
  file.remove(path)
  
}
