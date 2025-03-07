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
mkdirs('output/landsat_sensor_calibration/')

lsat.dt <- fread('data/lsat_extracts/lsat_arctic_clean_data_25k_sites_1999to2021_C2.csv')

# length(unique(lsat.dt$sample.id))
# sample.ids <- unique(lsat.dt$sample.id)
# subsample.ids <- sample(sample.ids, 100)
# lsat.dt <- lsat.dt[sample.id %in% subsample.ids]


# derive spectral indices
si.vec <- c('ndvi', 'evi2', 'ndwi', 'nbr','ndmi')
for (i in si.vec){
  lsat.dt <- lsat_calc_spectral_index(lsat.dt, si = i)  
}

# cross calibrate bands and spectral indices
band.vec <- c('blue','green','red','nir','swir1','swir2')
band.si.vec <- c(band.vec, si.vec)

i = 'ndmi'

for (i in band.si.vec){
  lsat.dt <- lsat_calibrate_poly(lsat.dt, 
                                 band.or.si = i, 
                                 overwrite.col = F, 
                                 write.output = T,
                                 outdir = 'output/landsat_sensor_calibration/')
  print(paste0('finished: ',i))
}


# combine all regression coefficients and evaluation statistics 
filenames <- list.files('output/landsat_sensor_calibration/', pattern = 'evaluation', full.names = T)
xcal.eval.dt <- do.call("rbind", lapply(filenames, fread))
xcal.eval.dt
fwrite(xcal.eval.dt, 'output/landsat_sensor_calibration/arctic_xcal_regression_evaluation.csv')
unlink(filenames)

filenames <- list.files('output/landsat_sensor_calibration/', pattern = 'coef', full.names = T)
xcal.coef.dt <- do.call("rbind", lapply(filenames, fread, fill=T))
xcal.coef.dt
fwrite(xcal.coef.dt, 'output/landsat_sensor_calibration/arctic_xcal_regression_coefficients.csv')
unlink(filenames)


# END SCRIPT ========================================================================================================

# highlat.xcal.eval.dt <- fread('output/arctic/50k/arctic_xcal_regression_evaluation.csv')
# arctic.xcal.eval.dt <- fread('output/arctic/arctic_xcal_regression_evaluation.csv')
# 
# xx <- highlat.xcal.eval.dt[arctic.xcal.eval.dt, on = c('satellite','band.or.si')]
# 
# plot(xx$r2, xx$i.r2, ylab='high lat', xlab='arctic')
# abline(a=0, b=1)
# cor.test(xx$r2, xx$i.r2)
# mean(xx$r2, na.rm = T)
# mean(xx$i.r2, na.rm = T)
# t.test(xx$r2, xx$i.r2, paired = T)
# 
# plot(xx$xcal.bias.pcnt, xx$i.xcal.bias.pcnt, ylab='high lat', xlab='arctic')
# abline(a=0, b=1)
# cor.test(xx$xcal.bias.pcnt, xx$i.xcal.bias.pcnt)
# t.test(xx$xcal.bias.pcnt, xx$i.xcal.bias.pcnt, paired = T)
# 
# plot(xx$uncal.bias.pcnt, xx$i.uncal.bias.pcnt, ylab='high lat', xlab='arctic')
# abline(a=0, b=1)
# cor.test(xx$uncal.bias.pcnt, xx$i.uncal.bias.pcnt)
# t.test(xx$uncal.bias.pcnt, xx$i.uncal.bias.pcnt, paired = T)
# 
# arctic.xcal.eval.dt[, .(r2.avg = round(mean(r2),2), 
#                         r2.sd = round(sd(r2),2),
#                         uncal.bias.pcnt.avg = round(mean(abs(uncal.bias.pcnt)),1), 
#                         uncal.bias.pcnt.sd = round(sd(abs(uncal.bias.pcnt)),1), 
#                         xcal.bias.pcnt.avg = round(mean(abs(xcal.bias.pcnt)),1), 
#                         xcal.bias.pcnt.sd = round(sd(abs(xcal.bias.pcnt)),1)), by = 'satellite']

