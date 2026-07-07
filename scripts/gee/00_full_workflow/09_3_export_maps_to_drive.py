"""

DESCRIPTION: Export final maps to Google Drive for external use

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:

TO-DO:

"""

import ee
from utils import params

try:
    ee.Initialize(project=params.ee_project)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=params.ee_project)

# =========================
# 1. SET UP ===============
# =========================

# 1.0 ----- PARAMETERS -----

# Model options
year = 2020

# Visualization options
biomass_max = 2000

# Reprojection options
crs = params.crs
scale = params.scale

# 1.1 ----- READ IN DATA -----

plant_biomass_percentiles = ee.Image("projects/arctic-biomass-mapping/assets/modeled_final/total_biomass_2020_j_index_v20240514")
woody_biomass_percentiles = ee.Image("projects/arctic-biomass-mapping/assets/modeled_final/woody_biomass_2020_j_index_v20240514")
water_snow_ice = ee.Image("projects/arctic-biomass-mapping/assets/modeled_final/water_snow_ice").rename('mask_val')
canopy_height = ee.ImageCollection("projects/meta-forest-monitoring-okw37/assets/CanopyHeight").mosaic()

# 1.2 ----- TREE MASK SET-UP -----

# Fill currently masked areas with low value to ensure we include Greeland, barrier islands etc.
canopy_height = canopy_height.unmask(0).setDefaultProjection('EPSG:4326')

# Create mask
canopy_height_mask = canopy_height.lt(5).rename('mask_val')

# =========================
# 2. PARTITION ============
# =========================

# 2.0 ----- MEDIAN -----

plant_biomass = plant_biomass_percentiles.select('total_biomass').rename('biomass_gm2')
woody_biomass = woody_biomass_percentiles.select('woody_biomass').rename('biomass_gm2')

# 2.1 ----- UNCERTAINTY BOUNDS -----

plant_biomass_lwr = plant_biomass_percentiles.select('total_biomass_lwr').rename('biomass_gm2')
plant_biomass_upr = plant_biomass_percentiles.select('total_biomass_upr').rename('biomass_gm2')
woody_biomass_lwr = woody_biomass_percentiles.select('woody_biomass_lwr').rename('biomass_gm2')
woody_biomass_upr = woody_biomass_percentiles.select('woody_biomass_upr').rename('biomass_gm2')

# =========================
# 3. WOODY DOMINANCE ======
# =========================

# 3.0 ----- CALCULATE WOODY DOMINANCE -----
# Note: where plant and woody plant biomass are both zero, woody plant dominance is also zero

woody_dominance = woody_biomass.divide(plant_biomass).multiply(100).rename('dominance_percent')

# =========================
# 4. FINALIZE MAPS ========
# =========================

# 4.0 ----- MASK WATER/SNOW/ICE -----

# Tidy water/snow/ice mask
water_snow_ice = water_snow_ice.multiply(0).clip(ee.Geometry.Polygon([-180, 81.86, 0, 81.86, 180, 81.86, 180, 45, 0, 45, -180, 45], None, False))

# Add water/snow/ice zeros
plant_biomass_lwr = ee.ImageCollection([plant_biomass_lwr, water_snow_ice.rename('biomass_gm2')]).mosaic()
plant_biomass = ee.ImageCollection([plant_biomass, water_snow_ice.rename('biomass_gm2')]).mosaic()
plant_biomass_upr = ee.ImageCollection([plant_biomass_upr, water_snow_ice.rename('biomass_gm2')]).mosaic()
woody_biomass_lwr = ee.ImageCollection([woody_biomass_lwr, water_snow_ice.rename('biomass_gm2')]).mosaic()
woody_biomass = ee.ImageCollection([woody_biomass, water_snow_ice.rename('biomass_gm2')]).mosaic()
woody_biomass_upr = ee.ImageCollection([woody_biomass_upr, water_snow_ice.rename('biomass_gm2')]).mosaic()
woody_dominance = ee.ImageCollection([woody_dominance, water_snow_ice.rename('dominance_percent')]).mosaic()

# 4.1 ----- TIDY CANOPY HEIGHT MASK -----

# Restrict canopy height mask to study area
canopy_height_mask = canopy_height_mask.updateMask(plant_biomass)

# =========================
# 5. EXPORT ===============
# =========================

task = ee.batch.Export.image.toDrive(
    image=plant_biomass_lwr.round().uint16().unmask(65535),  # Convert to int16 for easier storage
    description='plant_biomass_' + str(year) + '_p2_5',
    folder='final_maps',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    fileDimensions=89600,  # https://gis.stackexchange.com/questions/359974/calculating-max-shardsize-and-filedimensions-to-avoid-tiling-large-raster-export/360062#360062
    maxPixels=1e13,
    formatOptions={  # https://developers.google.com/earth-engine/guides/exporting_images
        'cloudOptimized': True,
        'noData': 65535
    }
)
task.start()

task = ee.batch.Export.image.toDrive(
    image=plant_biomass.round().uint16().unmask(65535),  # Convert to int16 for easier storage
    description='plant_biomass_' + str(year) + '_p50',
    folder='final_maps',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    fileDimensions=89600,  # https://gis.stackexchange.com/questions/359974/calculating-max-shardsize-and-filedimensions-to-avoid-tiling-large-raster-export/360062#360062
    maxPixels=1e13,
    formatOptions={  # https://developers.google.com/earth-engine/guides/exporting_images
        'cloudOptimized': True,
        'noData': 65535
    }
)
task.start()

task = ee.batch.Export.image.toDrive(
    image=plant_biomass_upr.round().uint16().unmask(65535),  # Convert to int16 for easier storage
    description='plant_biomass_' + str(year) + '_p97_5',
    folder='final_maps',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    fileDimensions=89600,  # https://gis.stackexchange.com/questions/359974/calculating-max-shardsize-and-filedimensions-to-avoid-tiling-large-raster-export/360062#360062
    maxPixels=1e13,
    formatOptions={  # https://developers.google.com/earth-engine/guides/exporting_images
        'cloudOptimized': True,
        'noData': 65535
    }
)
task.start()

task = ee.batch.Export.image.toDrive(
    image=woody_biomass_lwr.round().uint16().unmask(65535),  # Convert to int16 for easier storage
    description='woody_plant_biomass_' + str(year) + '_p2_5',
    folder='final_maps',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    fileDimensions=89600,  # https://gis.stackexchange.com/questions/359974/calculating-max-shardsize-and-filedimensions-to-avoid-tiling-large-raster-export/360062#360062
    maxPixels=1e13,
    formatOptions={  # https://developers.google.com/earth-engine/guides/exporting_images
        'cloudOptimized': True,
        'noData': 65535
    }
)
task.start()

task = ee.batch.Export.image.toDrive(
    image=woody_biomass.round().uint16().unmask(65535),  # Convert to int16 for easier storage
    description='woody_plant_biomass_' + str(year) + '_p50',
    folder='final_maps',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    fileDimensions=89600,  # https://gis.stackexchange.com/questions/359974/calculating-max-shardsize-and-filedimensions-to-avoid-tiling-large-raster-export/360062#360062
    maxPixels=1e13,
    formatOptions={  # https://developers.google.com/earth-engine/guides/exporting_images
        'cloudOptimized': True,
        'noData': 65535
    }
)
task.start()

task = ee.batch.Export.image.toDrive(
    image=woody_biomass_upr.round().uint16().unmask(65535),  # Convert to int16 for easier storage
    description='woody_plant_biomass_' + str(year) + '_p97_5',
    folder='final_maps',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    fileDimensions=89600,  # https://gis.stackexchange.com/questions/359974/calculating-max-shardsize-and-filedimensions-to-avoid-tiling-large-raster-export/360062#360062
    maxPixels=1e13,
    formatOptions={  # https://developers.google.com/earth-engine/guides/exporting_images
        'cloudOptimized': True,
        'noData': 65535
    }
)
task.start()

task = ee.batch.Export.image.toDrive(
    image=woody_dominance.round().byte().unmask(255),  # Convert to uint8 for easier storage
    description='woody_plant_dominance_' + str(year),
    folder='final_maps',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    fileDimensions=89600,  # https://gis.stackexchange.com/questions/359974/calculating-max-shardsize-and-filedimensions-to-avoid-tiling-large-raster-export/360062#360062
    maxPixels=1e13,
    formatOptions={  # https://developers.google.com/earth-engine/guides/exporting_images
        'cloudOptimized': True,
        'noData': 255
    }
)
task.start()

task = ee.batch.Export.image.toDrive(
    image=canopy_height_mask.round().byte().unmask(255),  # Convert to uint8 for easier storage
    description='canopy_height_mask_lt5m',
    folder='final_maps',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    fileDimensions=89600,  # https://gis.stackexchange.com/questions/359974/calculating-max-shardsize-and-filedimensions-to-avoid-tiling-large-raster-export/360062#360062
    maxPixels=1e13,
    formatOptions={  # https://developers.google.com/earth-engine/guides/exporting_images
        'cloudOptimized': True,
        'noData': 255
    }
)
task.start()

raise Exception('stop')