"""

DESCRIPTION: Calculate median DOY for each season across Pan Arctic - this provides, for each season, a reasonable date to use to extract a CCDC segment

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:
- This script uses the original Arctic polygon without updated Greenland boundaries
- This is so the entire Greenland icesheet is not used in determining median season DOYs

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

# CCDC input options
ccdc_start_year = 1984
ccdc_end_year = 2023

# Output options
version = 'v20240207'  # Version identifier

# Reprojection options
crs = params.crs
scale = params.scale

# 1.1 ----- READ IN DATA -----

arctic = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/arctic_oroarctic_polygon')
seasonal_doys = ee.ImageCollection('projects/arctic-biomass-mapping/assets/seasonal_percentile_doys/seasonal_doys_tiles_' + str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + version)
seasonal_doy = seasonal_doys.mosaic()

# =========================
# 2. ANALYSIS =============
# =========================

doy_median = seasonal_doy.reduceRegions(
    reducer=ee.Reducer.median(),
    scale=scale,
    crs=crs,
    collection=arctic
)

# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.table.toAsset(
    collection=doy_median,
    description='seasonal_doys_median',
    assetId='projects/arctic-biomass-mapping/assets/seasonal_percentile_doys/seasonal_doys_median_' + str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + version
)
task.start()

raise Exception('stop')
