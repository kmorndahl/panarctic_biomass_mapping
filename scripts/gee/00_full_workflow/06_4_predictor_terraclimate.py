"""

DESCRIPTION: Create climate predictors

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

# 1.0 ----- READ IN DATA -----

# Read in TerraClimate data
# https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_TERRACLIMATE
terra_climate = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")

# 1.1 ----- PARAMETERS -----

# TerraClimate date range
ccdc_start_year = 1995
ccdc_end_year = 2022

# Projection (derived from input data)
proj = terra_climate.first().projection().getInfo()
crs = proj['wkt']
transform = proj['transform']

# =========================
# 2. ANALYSIS =============
# =========================

# 2.0 ----- CALCULATE VARIABLES -----

# Get number of years
climate_nyears = ee.Number(ccdc_end_year).subtract(ee.Number(ccdc_start_year)).add(1)

# Calculate average yearly temperature
def calc_tmean(i):
    return i.expression("tmmean = (b('tmmx') + b('tmmn'))/2.0").divide(10)  # Scale factor

climate_tmean = (terra_climate
    .filter(ee.Filter.calendarRange(ccdc_start_year, ccdc_end_year, "year"))
    .map(calc_tmean))

# Calculate summer warmth index
climate_swi = (climate_tmean
    .map(lambda i: i.multiply(i.gt(0)))
    .reduce(ee.Reducer.sum())
    .divide(climate_nyears)
    .rename('swi'))

# Calculate minimum January temperature
climate_tmin = (terra_climate
    .filter(ee.Filter.calendarRange(1, 1, "month"))
    .filter(ee.Filter.calendarRange(ccdc_start_year, ccdc_end_year, "year"))
    .select('tmmn')
    .reduce(ee.Reducer.mean())
    .divide(10)  # Scale factor
    .rename('tmin'))

# Calculate total annual precipitation
climate_tapmm = (terra_climate
    .filter(ee.Filter.calendarRange(ccdc_start_year, ccdc_end_year, "year"))
    .select('pr')
    .reduce(ee.Reducer.sum())
    .divide(climate_nyears)
    .rename('tapmm'))

# Combine
climate = (climate_swi
    .addBands(climate_tmin)
    .addBands(climate_tapmm))

# 2.1 ----- FILL NAs -----

terraclimate = ee.Image("projects/arctic-biomass-mapping/assets/predictors/terraclimate")
terraclimate_fill_gaps = terraclimate.focalMean(
    radius=terraclimate.projection().nominalScale().multiply(3).divide(2),
    kernelType='square',
    units='meters',
    iterations=1
)  # Fill some gaps to avoid NAs in cal/val data
terraclimate = terraclimate_fill_gaps.blend(terraclimate)
terraclimate = terraclimate.regexpRename('^', 'climate_4km_')

# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.image.toAsset(
    image=terraclimate,
    description='terraclimate_filled',
    assetId='projects/arctic-biomass-mapping/assets/predictors/terraclimate_filled',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    crs=crs,
    crsTransform=transform,
    maxPixels=1e12
)
task.start()

raise Exception('stop')