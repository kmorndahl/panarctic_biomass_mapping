# DESCRIPTION =============================================================================\
# This R script generates Landsat cross-calibration regression models for the Arctic by 
# extracting data for 10,000 random sample locations
# Author: Logan Berner 
# Institution: Northern Arizona University, School of Informatics, Computing, and Cyber Systems
# Date: 2022-03-27
# URL: https://github.com/logan-berner/lsatTS
# =========================================================================================

# Clear workspace
rm(list=ls())

# Load required R packages
require(sf)
require(geosphere)
require(LandsatTS)
require(rgee)
require(dplyr)
require(R.utils)
require(data.table)

sf::sf_use_s2(FALSE)
mkdirs('data/lsat_extracts/raw')

arctic.sf <- st_read('data/arctic_oroarctic_laea.shp')
plot(arctic.sf)


# CREATE SAMPLE POINTS =========================================================================================

# sample high-latitude study area
sample.pts <- st_sample(arctic.sf, size = 15000) %>% st_cast('POINT')
sample.pts <- sample.pts %>% st_sf %>% mutate(sample_id = paste0('pt', 10001:(10000+nrow(st_coordinates(sample.pts)))))
plot(sample.pts, col = 'blue')

# cluster sample points to speed data extraction from GEE
sample.pts.latlon <- st_transform(sample.pts, crs = st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
coords <- st_coordinates(sample.pts.latlon)
dist <- distm(coords) # distance matrix
clusters <- hclust(as.dist(dist), method="complete")
sample.pts$cluster <- cutree(clusters, h=1000000)
length(unique(sample.pts$cluster)) # number of clusters

# quick map 
plot(sample.pts)

# write out shapefile
st_write(sample.pts, 'data/lsat_extracts/lsat_arctic_random_sample_pts_n15000.shp')

# EXTRACT LANDSAT DATA FOR SAMPLE POINTS =========================================================================
ee_Initialize()
task_list <- lsat_export_ts(pixel_coords_sf = sample.pts, chunks_from = 'cluster', 
                            start_doy = 121, end_doy = 273, 
                            start_date = '1999-05-01', end_date = '2022-09-30',
                            file_prefix = 'arctic', drive_export_dir = 'earth_engine')

# END SCRIPT ========================================================================================================