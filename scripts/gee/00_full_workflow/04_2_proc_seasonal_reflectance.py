"""

DESCRIPTION: Produce seasonal composites from CCDC models based seasonal DOY images

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:

- Reprojection of individual thematic predictor sets NOT necessary
- CCDC memory use discussion: https://groups.google.com/g/google-earth-engine-developers/c/shIVv-ZcIzU/m/oCMXaADLBQAJ
- segmentFindStrategy
  - For clipping_method = tiles
    - Choose 'next' -- in the case where a disturbance occurs during the snow-free season 'next' allows the disturbance to be detected in the year it occurs
  - For clipping_method = calval
    - Choose 'closest'
    - For field sites from the beginning or end of the time series e.g., 2022, there is no 'next' segment so this allows 'previous' to be selected
    - In other instances, we will discard plots if they occur during a break between segments
    - These are likely due to disturbance and might mess up model fits
    - We will make exceptions for harvest dates just before the start or just after the end of the time series
    - In these cases, the lack of segment might be more to do with data availability than disturbance

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
year = 2020  # Select year for predictions

# CCDC input options
ccdc_start_year = 1984
ccdc_end_year = 2023

# CCDC processing options
extrapolate_max_days = 120  # Number of days to extrapolate beyond the start and end of a CCDC segment, helps fill in gaps before the first segment, after the last segment, and between segments

# Output data options
crs = params.crs
scale = params.scale
version = 'v20240207'  # Input version identifier
clipping_method = 'tiles'  # Choose 'tiles' or 'calval'
out_path = 'projects/arctic-biomass-mapping/assets/seasonal_modeled_reflectance/'

# Strategy to use for selecting CCDC segment when there is no segment for the specified date
# Choose 'previous', 'next' or 'closest'
segment_find_strategy = 'closest' if clipping_method == 'calval' else 'next'

# 1.1 ----- EXISTING MODULES/FUNCTIONS -----

from ccdc_tools import dates as utils_dates

from utils import misc as fun_misc
from utils import refl as fun_refl
from utils import temporal_segmentation as temporal_seg
Segments = temporal_seg.Segments

# FUNCTION: convert_dates
# USE: Convert from DOY to ee.Date using utils_dates.convertDate
# NOTES: Mapped over dictionary (applied via ee.Dictionary.map in Python)
# AUTHOR: Katie Orndahl
# LAST UPDATE: 11-10-2024

def convert_dates(key, value):
    frac_date = ee.Number(year).add(ee.Number(value).divide(365.25))  # Convert from DOY to fractional year
    return utils_dates.convertDate({'inputFormat': 1, 'inputDate': frac_date, 'outputFormat': 4})

# 1.2 ----- READ IN DATA -----

seasonal_doys = ee.ImageCollection('projects/arctic-biomass-mapping/assets/seasonal_percentile_doys/seasonal_doys_tiles_' + str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + version)
ccdc_fits = ee.ImageCollection('projects/arctic-biomass-mapping/assets/CCDC_tiles/CCDC_C2_SR_tiles_' + str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + version)
doy_median = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/seasonal_percentile_doys/seasonal_doys_median_' + str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + version)

# =========================
# 2. ANALYSIS =============
# =========================

# 2.0 ----- MOSAIC DATA -----

ccdc_fit = ccdc_fits.mosaic()
seasonal_doy = seasonal_doys.mosaic()

# 2.1 ----- MODEL REFLECTANCE -----

# Get median seasonal DOYs
# This provides, for each season, a reasonable date to use to extract a CCDC segment
doy_median = doy_median.first().toDictionary(['doy_earlySummer', 'doy_endSnowfree', 'doy_lateSummer', 'doy_peakSummer', 'doy_startSnowfree'])

# Convert DOYs to fractional dates
doy_median_frac = doy_median.map(convert_dates)

# Get segments
segments = Segments(ccdc_fit, date_format=1)  # 1 = provide dates as fractional years

# Model CCDC reflectance and rename bands with seasonal suffix
# Segments are selected based on average DOY for the given season
start_snowfree_fit = fun_refl.fitT(
    seasonal_doy.select('doy_startSnowfree'),
    year,
    segments.find_by_date(ee.Date(doy_median_frac.get('doy_startSnowfree')), strategy=segment_find_strategy),
    extrapolate_max_days
).regexpRename('$', '_startSnowfree').regexpRename('^', 'spectral_')

early_summer_fit = fun_refl.fitT(
    seasonal_doy.select('doy_earlySummer'),
    year,
    segments.find_by_date(ee.Date(doy_median_frac.get('doy_earlySummer')), strategy=segment_find_strategy),
    extrapolate_max_days
).regexpRename('$', '_earlySummer').regexpRename('^', 'spectral_')

peak_summer_fit = fun_refl.fitT(
    seasonal_doy.select('doy_peakSummer'),
    year,
    segments.find_by_date(ee.Date(doy_median_frac.get('doy_peakSummer')), strategy=segment_find_strategy),
    extrapolate_max_days
).regexpRename('$', '_peakSummer').regexpRename('^', 'spectral_')

late_summer_fit = fun_refl.fitT(
    seasonal_doy.select('doy_lateSummer'),
    year,
    segments.find_by_date(ee.Date(doy_median_frac.get('doy_lateSummer')), strategy=segment_find_strategy),
    extrapolate_max_days
).regexpRename('$', '_lateSummer').regexpRename('^', 'spectral_')

end_snowfree_fit = fun_refl.fitT(
    seasonal_doy.select('doy_endSnowfree'),
    year,
    segments.find_by_date(ee.Date(doy_median_frac.get('doy_endSnowfree')), strategy=segment_find_strategy),
    extrapolate_max_days
).regexpRename('$', '_endSnowfree').regexpRename('^', 'spectral_')

# Combine modeled reflectance
modeled_reflectance = start_snowfree_fit.addBands(early_summer_fit) \
                                        .addBands(peak_summer_fit) \
                                        .addBands(late_summer_fit) \
                                        .addBands(end_snowfree_fit) \
                                        .int16()

print('Seasonal modeled reflectance', modeled_reflectance.getInfo())

# =========================
# 3. EXPORT ===============
# =========================

# Ensure output ImageCollection exists
fun_misc.createCollectionIfNotExists(out_path.rstrip('/'))

task = ee.batch.Export.image.toAsset(
    image=modeled_reflectance,
    description='seasonal_modeled_reflectance_' + str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + version + '_' + str(year),
    assetId=out_path + 'seasonal_modeled_reflectance_' + str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + version + '_' + str(year),
    scale=scale,
    crs=crs,
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    maxPixels=1e13
)
task.start()

raise Exception('stop')
