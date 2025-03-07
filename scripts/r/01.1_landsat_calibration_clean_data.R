# DESCRIPTION =============================================================================\
# This R script generates Landsat cross-calibration regression models for the Arctic by 
# extracting data for 25,000 random sample locations
# Author: Logan Berner 
# Institution: Northern Arizona University, School of Informatics, Computing, and Cyber Systems
# Date: 2023-02-06
# =========================================================================================
# Clear workspace
rm(list=ls())
setwd('data/lsat_extracts/')

# Load required R packages
require(LandsatTS)
require(data.table)
require(R.utils)

sf::sf_use_s2(FALSE)
mkdirs('clean/')

data.files <- list.files('raw/', full.names = T)
n.data.files <- length(data.files)

for (i in 1:n.data.files){
  lsat.dt <- fread(data.files[i])
  lsat.dt <- lsat_format_data(lsat.dt)
  lsat.dt <- lsat_clean_data(lsat.dt)
  fwrite(lsat.dt, file = paste0('clean/',i,'.csv'))
  print(i/n.data.files)
}

# Combine cleaned files and write out as one CSV 
# data.files <- list.files('data/arctic_sample_25000_points/clean/', full.names = T)
# lsat.dt <- do.call("rbind", lapply(data.files, fread))
lsat.dt <- fwrite(lsat.dt, 'lsat_arctic_clean_data_25k_sites_1999to2021_C2.csv')

# END SCRIPT ========================================================================================================