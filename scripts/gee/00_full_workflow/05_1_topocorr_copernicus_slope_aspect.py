"""

DESCRIPTION: Derive slope and aspect from Copernicus DEM -- saves compute in later scripts

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

# Reprojection options
copernicus_crs = 'EPSG:3857'  # Native CRS of Copernicus DEM
scale = params.scale

# Output options
no_data_val = -9999

# 1.1 ----- READ IN DATA -----

# Copernicus DEM, 30 m resolution
dem_ic = ee.ImageCollection('COPERNICUS/DEM/GLO30');

# Need to reproject DEM to calculate slope:
# https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_DEM_GLO30
# https://groups.google.com/g/google-earth-engine-developers/c/KwKvcS_3QzI/m/pqLdorsqCQAJ
dem = dem_ic.mosaic() \
        .setDefaultProjection(copernicus_crs, None, scale) \
        .select('DEM') \
        .rename('elevation')

# ROI
arctic_poly = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/arctic_oroarctic_coast_buffer_laea_final')

# =========================
# 2. ANALYSIS =============
# =========================

# 2.0 ----- SMOOTH DEM -----

# Smooth the DEM using a low-pass kernel
# reduceNeighborhood is more efficient, see here: https://groups.google.com/g/google-earth-engine-developers/c/sfvZlGQg5yw/m/lgvKV0b0BQAJ
boxcar = ee.Kernel.square(radius=3, units='pixels', normalize=True)
dem_s = ee.Image(dem).reduceNeighborhood(ee.Reducer.mean(), boxcar)  # Confirmed this produces identical results to .convolve()

# 2.1 ----- CREATE TERRAIN LAYERS -----

slp = ee.Terrain.slope(dem_s) \
                .rename('slope') \
                .clip(arctic_poly)

asp = ee.Terrain.aspect(dem_s) \
                .rename('aspect') \
                .clip(arctic_poly)

print('Slope', slp.getInfo())
print('Aspect', asp.getInfo())

# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.image.toCloudStorage(
    image=slp,
    description='slp',
    bucket='biomass_predictors',
    fileNamePrefix='topographic/copernicus_30m/copernicus_30m_slope/copernicus_slp_',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False), 
    scale=scale,
    crs=copernicus_crs,
    maxPixels=1e13,
    formatOptions={
        'cloudOptimized': True,
        'noData': no_data_val
    }
)
task.start()

task = ee.batch.Export.image.toCloudStorage(
    image=asp,
    description='asp',
    bucket='biomass_predictors',
    fileNamePrefix='topographic/copernicus_30m/copernicus_30m_aspect/copernicus_asp_',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False), 
    scale=scale,
    crs=copernicus_crs,
    maxPixels=1e13,
    formatOptions={
        'cloudOptimized': True,
        'noData': no_data_val
    }
)
task.start()

raise Exception('stop')
