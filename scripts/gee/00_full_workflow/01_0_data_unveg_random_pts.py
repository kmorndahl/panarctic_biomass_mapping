"""

DESCRIPTION: Generate zero biomass observations for future modeling

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:
- This script uses the original Arctic polygon without updated Greenland boundaries
- This is so absences are more spread out, and not concentrated on the Greenland icesheet

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

# 1.0 ----- READ IN DATA -----

wc = ee.ImageCollection('ESA/WorldCover/v200').first()  # WorldCover dataset
arctic_poly = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/arctic_oroarctic_polygon')

# 1.1 ----- PARAMETERS -----

# Reprojection options
crs = params.crs
scale = params.scale

# =========================
# 2. ANALYSIS =============
# =========================

# 2.0 ----- PREPARE LAND COVER IMAGE -----

# Restrict to barren type land cover types only
built = wc.eq(50).selfMask()
barrens = wc.eq(60).selfMask()
snow_ice = wc.eq(70).selfMask()
water = wc.eq(80).selfMask()
lc = ee.Image(ee.ImageCollection([built, barrens, snow_ice, water]).mosaic().multiply(wc))

# Clip to arctic polygon
lc = lc.clip(arctic_poly)

# Visualize

# 2.1 ----- EXTRACT SAMPLES -----

random_pts = lc.stratifiedSample(
    numPoints=500,
    seed=1908,
    projection=crs,
    scale=scale,
    region=ee.Geometry.BBox(-180, 45, 179.9999, 89.9999),
    geometries=True
)

# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.table.toDrive(
    collection=random_pts,
    description='barren_built_snow_ice_water_random_pts_scale30m',
    folder='land_cover_random_pts',
    fileFormat='SHP'
)
task.start()

task = ee.batch.Export.table.toAsset(
    collection=random_pts,
    description='barren_built_snow_ice_water_random_pts_scale30m',
    assetId='projects/arctic-biomass-mapping/assets/field_data/barren_built_snow_ice_water_random_pts_scale30m'
)
task.start()

raise Exception('stop')

