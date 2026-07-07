"""

DESCRIPTION: Create auxillary modeled seasonal reflectance predictors

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:
- segmentFindStrategy
  - For clipping_method = tiles
    - Choose 'next' -- in the case where a disturbance occurs during the snow-free season 'next' allows the disturbance to be detected in the year it occurs
  - For clipping_method = calval
    - Choose 'closest'
    - For field sites from the beginning or end of the time series e.g., 2022, there is no 'next' segment so this allows 'previous' to be selected
    - In other instances, we will discard plots if they occur during a break between segments
    - These are likely due to disturbance and might mess up model fits
    - We will make exceptions for harvest dates just before the start or just after the end of the time series
    - In these cases, the lack of segment might be more to due with data availability than disturbance

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

# 1.0 ----- PARAMETERS ----

# CCDC date range
start_year = 1984
end_year = 2023

# Model options
year = 2000  # Select year for predictions
date = ee.Date.fromYMD(year, 7, 31)  # Create date for extracting segment, choose mid-summer-ish for month and day

# CCDC processing options
extrapolateMaxDays = 120  # Number of days to extrapolate beyond the start and end of a CCDC segment, helps fill in gaps before the first segment, after the last segment, and between segments

# Output data options
ccdc_version = 'v20240207'
crs = params.crs
scale = params.scale  # Scale for reprojecting predictors
clipping_method = 'tiles'  # Choose 'tiles' or 'calval'

# Texture analysis options
texture_bands = ['spectral_NDVI_peakSummer']  # Chosen from review of literature
texture_metrics = ['.*_var']  # Chosen from review of literature
texture_radius = 1  # Chosen from review of literature -- 3x3 window most common

# Strategy to use for selecting CCDC segment when there is no segment for the specified date
# Choose 'previous', 'next' or 'closest'
segmentFindStrategy = 'closest' if clipping_method == 'calval' else 'next'

# 1.1 ----- READ IN DATA ----

ccdc_fit = ee.ImageCollection(
    'projects/arctic-biomass-mapping/assets/CCDC_tiles/CCDC_C2_SR_tiles_' + str(start_year) + '_' + str(end_year) + '_' + ccdc_version
).mosaic()

refl = ee.Image(
    'projects/arctic-biomass-mapping/assets/seasonal_modeled_reflectance/seasonal_modeled_reflectance_' + str(start_year) + '_' + str(end_year) + '_' + ccdc_version + '_' + str(year)
)

refl_tc = ee.ImageCollection(
    'projects/arctic-biomass-mapping/assets/seasonal_modeled_reflectance/seasonal_modeled_reflectance_' + str(start_year) + '_' + str(end_year) + '_' + ccdc_version + '_' + str(year) + '_topographic_correction'
)

# 1.2 ----- TIDY DATA ----

refl_tc = refl_tc.mosaic()  # Mosaic topographically corrected collection
refl = refl.select('.*_ND.*|.*NBR.*')  # Get normalized indices from original image
modeled_reflectance = refl_tc.addBands(refl)  # Combine

# 1.3 ----- MODULES/FUNCTIONS ----

# Load existing modules
from utils import temporal_segmentation as temporal_seg
Segments = temporal_seg.Segments

# Load project function modules (working directory must be scripts/gee/)
from utils import refl as fun_refl

# =========================
# 2. ANALYSIS =============
# =========================

# 2.0 ----- SET UP SEGMENT ----

# Get segments
segments = Segments(ccdc_fit, date_format=1)

# Get segment that matches prediction date
segment = segments.find_by_date(date, strategy=segmentFindStrategy)

# 2.1 ----- TEXTURE ----

# Calculate texture
texture = fun_refl.calculateTexture(modeled_reflectance, texture_radius, texture_bands, texture_metrics, crs, scale)

# Fill gaps
texture_fill_gaps = texture.focalMean(
    radius=scale * 3 // 2,
    kernelType='square',
    units='meters',
    iterations=1
)  # Fill some gaps to avoid NAs in cal/val data
texture = (texture_fill_gaps.blend(texture)  # Fill some gaps to avoid NAs in cal/val data
                            .uint16())  # Set data type

# Report
print('NDVI texture', texture.getInfo())

# 2.2 ----- SLOPE ----

# Calculate NDVI slope
slope = fun_refl.NDVIslope(segment)

# Fill gaps
slope_fill_gaps = slope.focalMean(
    radius=scale * 3 // 2,
    kernelType='square',
    units='meters',
    iterations=1
)  # Fill some gaps to avoid NAs in cal/val data
slope = (slope_fill_gaps.blend(slope)  # Fill some gaps to avoid NAs in cal/val data
                        .int32()  # Set data type
                        .setDefaultProjection(crs=crs, scale=scale))  # Set default projection

# Report
print('NDVI slope', slope.getInfo())

# =========================
# 3. EXPORT ===============
# =========================

# Texture
task = ee.batch.Export.image.toAsset(
    image=texture,
    description='seasonal_modeled_reflectance_NDVI_texture_' + str(year) + '_' + ccdc_version,
    assetId='projects/arctic-biomass-mapping/assets/predictors/seasonal_modeled_reflectance_NDVI_texture_' + str(year) + '_' + ccdc_version,
    pyramidingPolicy='mean',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    crs=crs,
    scale=scale,
    maxPixels=1e13
)
task.start()

# Slope
task_slope = ee.batch.Export.image.toAsset(
    image=slope,
    description='seasonal_modeled_reflectance_NDVI_slope_' + str(year) + '_' + ccdc_version,
    assetId='projects/arctic-biomass-mapping/assets/predictors/seasonal_modeled_reflectance_NDVI_slope_' + str(year) + '_' + ccdc_version,
    pyramidingPolicy='mean',
    region=ee.Geometry.Polygon([-180, 83, 0, 83, 180, 83, 180, 45, 0, 45, -180, 45], None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    maxPixels=1e13
)
task_slope.start()

raise Exception('stop')