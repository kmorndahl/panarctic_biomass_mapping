"""

DESCRIPTION: Export example areas for hillshade graphic

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
# 1. SET-UP ===============
# =========================

# 1.0 ----- PARAMETERS -----

roi_name = 'ural'  # Choose 'ns' or 'ural'

# Reprojection options
scale = params.scale

# 1.1 ----- READ IN DATA -----

plant_biomass_percentiles = ee.Image("projects/arctic-biomass-mapping/assets/modeled_final/total_biomass_2020_j_index_v20240514")
water_snow_ice = ee.Image("projects/arctic-biomass-mapping/assets/modeled_final/water_snow_ice").rename('mask_val')
dem = ee.Image("UMN/PGC/ArcticDEM/V3/2m_mosaic")

ural = ee.Geometry.Polygon(
    [[[66.7356702920265, 68.22969954622916],
      [66.7356702920265, 68.12555897974919],
      [67.01856824124525, 68.12555897974919],
      [67.01856824124525, 68.22969954622916]]], None, False)

ns = ee.Geometry.Polygon(
    [[[-141.51688518059237, 69.4947680051644],
      [-141.51688518059237, 69.44106329224476],
      [-141.34179057610018, 69.44106329224476],
      [-141.34179057610018, 69.4947680051644]]], None, False)

roi = ee.Algorithms.If(ee.String(roi_name).equals('ns'),
                       ns,
                       ural)

# =========================
# 2. TIDY DATA ============
# =========================

# 2.0 ----- PARTITION DATA -----

# Get median
plant_biomass = plant_biomass_percentiles.select('total_biomass').rename('plant_biomass_gm2')

# 2.1 ----- ADD ZEROS -----

water_snow_ice = water_snow_ice.multiply(0).clip(ee.Geometry.Polygon([-180, 81.86, 0, 81.86, 180, 81.86, 180, 45, 0, 45, -180, 45], None, False))
plant_biomass = ee.ImageCollection([plant_biomass, water_snow_ice.rename('plant_biomass_gm2')]).mosaic()

# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.image.toDrive(
    image=plant_biomass,
    description='hillshade_' + roi_name + '_biomass',
    folder='hillshade',
    region=roi,
    scale=scale,
    crs='EPSG:3857'
)
task.start()

task = ee.batch.Export.image.toDrive(
    image=dem,
    description='hillshade_' + roi_name + '_arcticdem',
    folder='hillshade',
    region=roi,
    scale=scale,
    crs='EPSG:3857'
)
task.start()

raise Exception('stop')