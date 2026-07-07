"""

DESCRIPTION: Create transects along topographic gradients and extract biomass

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

# Region options
roi_name = 'ns'  # Choose 'ural' or 'ns'

# Visualization options
biomass_max = 2000

# Reprojection options
crs = params.crs
scale = params.scale

# 1.1 ----- READ IN DATA -----

ural = ee.Geometry.LineString(
    [[66.85925157289955, 68.16970307308829],
     [66.88740403871986, 68.18680415314952]])

ns = ee.Geometry.LineString(
    [[-141.44121332087335, 69.45968016053793],
     [-141.42130060114678, 69.46214945740564]])

roi = ee.Algorithms.If(ee.String(roi_name).equals('ural'), ural, ns)
plant_biomass_percentiles = ee.Image("projects/arctic-biomass-mapping/assets/modeled_final/total_biomass_2020_j_index_v20240514")
water_snow_ice = ee.Image("projects/arctic-biomass-mapping/assets/modeled_final/water_snow_ice").rename('mask_val')

# Copernicus DEM, 30 m resolution
dem_tc = (ee.ImageCollection("COPERNICUS/DEM/GLO30")
          .mosaic()
          .setDefaultProjection('EPSG:3857', None, 30)  # Explanation on setting default projection here https://twitter.com/jstnbraaten/status/1494038930643042309
          .select('DEM')
          .rename('elevation'))

# 1.2 ----- TIDY DATA -----

# Tidy biomass data
plant_biomass = plant_biomass_percentiles.select('total_biomass').rename('plant_biomass_gm2')
plant_biomass = ee.ImageCollection([plant_biomass, water_snow_ice.rename('plant_biomass_gm2')]).mosaic()

# Combine biomass and elevation data
biomass_elev = plant_biomass.addBands(dem_tc)

# =========================
# 2. ANALYSIS =============
# =========================

# 2.0 ----- CREATE TRANSECT -----

# Add points
pts = ee.FeatureCollection.randomPoints(roi, 1000, 65)

# 2.1 ----- EXTRACT DATA -----

biomass_elev_transect_fc = biomass_elev.reduceRegions(
    collection=pts,
    reducer=ee.Reducer.first(),
    scale=scale,
    crs=crs
)

# =========================
# 3. EXPORT ===============
# =========================

# points .csv
task = ee.batch.Export.table.toDrive(
    collection=biomass_elev_transect_fc,
    description='biomass_elev_transect_fc_' + roi_name,
    folder='biomass_topography_transect'
)
task.start()

# points .shp
task = ee.batch.Export.table.toDrive(
    collection=biomass_elev_transect_fc,
    description='biomass_elev_transect_fc_' + roi_name,
    folder='biomass_topography_transect',
    fileFormat='SHP'
)
task.start()

# line .shp
task = ee.batch.Export.table.toDrive(
    collection=roi,
    description='line_' + roi_name,
    folder='biomass_topography_transect',
    fileFormat='SHP'
)
task.start()

raise Exception('stop')