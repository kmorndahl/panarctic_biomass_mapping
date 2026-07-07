"""

DESCRIPTION: Apply CCDC models to produce seasonal composites

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:
- Start and end of the snowfree season DOY composites tend to be blocky due to the DOY restrictions associated with the MGRS tiles
  - However, these restrictions are necessary to get good CCDC fits
  - It does not appear that this blocky-ness impacts the modeled reflectance
  - More from Matt: "Remember that spring and fall dates are further clipped by the absolute DOY filter, which is applied based on latitude.
    It excluded dates when the peak sun elevation is less than 40 degrees in spring, and less than 25 degrees in fall.
    Mainly this is to exclude some weird garbage observations that would otherwise creep in, especially but not entirely on north-facing shadowy slopes.
    These garbage observations are mostly undetected snow and/or deep shadows and
    including them in the CCDC fits tends to yield some noisy results for the spring and fall windows especially.
    Using the peak sun angle is a way to come up with a general rule that varies with latitude.
    But, I have seen some evidence (fairly lush 'spring' composites) that these thresholds are not working in some southern coastal areas.
    It is also common that it has not snowed yet by the time the fall window closes.
    But at that point the sun angle is so low that the shadows are long and signal to noise ratio is low."
- Peak summer and other phenology DOY resources for comparison:
  - https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.14638
  - https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020JG006094
  - https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2486.2007.01529.x
  - https://tc.copernicus.org/articles/12/3373/2018/
  - https://iopscience.iop.org/article/10.1088/1748-9326/6/3/035502/pdf
  - https://www.sciencedirect.com/science/article/pii/S0924271618303010
  - https://www.cambridge.org/core/journals/annals-of-glaciology/article/vegetation-phenology-in-greenland-and-links-to-cryospheric-change/E5E1F5CE72E7A228B7C35B37D7B140C3
  - https://link.springer.com/article/10.1186/1472-6785-7-9#Sec13
- Finding an appropriate number of years/days to iterate over required some trial and error
  - Settled on 8 years, 77 days, 616 total model fits
  - Out of memory issues crop up when using over 700 model fits
- Choosing 'previous' as the segmentFindStrategy prevents peak summer identification from being affected by disturbance
- Index to use for determining peak summer greenness
  - In general, EVI is probably better suited for this
  - However, EVI is very prone to topographic illumination artifacts
  - Because we are applying the topographic correction post-hoc, we need to use NDVI here
  - Normalized differenced indices do not need to be topographically corrected
  - References:
    - https://www.sciencedirect.com/science/article/pii/S016819231830087X?casa_token=3RWQXv_WbtIAAAAA:YiK-weXtaIbQ-IADIdm0kqgJiixJXlRClT9D8PIQISA2ihuub3p6vzRLQsEIR3qmAmvGZ_z0jf8
    - https://www.mdpi.com/2673-7418/3/1/12
    - https://daac.ornl.gov/ABOVE/guides/Annual_Seasonality_greenness.html
    - https://www.sciencedirect.com/science/article/pii/S0303243419312164
    - https://www.sciencedirect.com/science/article/pii/S0034425714001011
    - https://www.mdpi.com/2072-4292/12/14/2290
- CCDC model fit band explanations
 * tStart: The start date of each model segment.
 * tEnd: The end date of each model segment.
 * tBreak: The model break date if a change is detected.
 * numObs: The number of observations used in each model segment.
 * changeProb: A numeric value representing the multi-band change probability.
 * _coefs: The regression coefficients for each of the bands in the image collection.
 * _rmse: The model root-mean-square error for each segment and input band.
 * _magnitude: For segments with changes detected, this represents the normalized residuals during the change period.

TO-DO:

"""

import ee
from utils import params
import sys

try:
    ee.Initialize(project=params.ee_project)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=params.ee_project)

# =========================
# 1. SET UP ===============
# =========================

# 1.0 ----- PARAMETERS -----

# CCDC processing options
segment_find_strategy = 'previous'  # Strategy to use for selecting CCDC segment when there is no segment for the specified date
peak_summer_index = 'NDVI'  # Index to use for determining peak summer greenness

# CCDC date range
ccdc_start_year = 1984
ccdc_end_year = 2023
start_MD = '06-15'
end_MD = '08-31'
map_year_start = ee.Number(ccdc_start_year).divide(5).ceil().multiply(5)
map_year_end = ee.Number(ccdc_end_year).divide(5).floor().multiply(5)
map_years = ee.List.sequence(map_year_start, map_year_end, 5)

# DOY band names
start_snowfree_name = 'doy_p003'
end_snowfree_name = 'doy_p097'
extrapolate_max_days = 120  # Number of days to extrapolate beyond the start and end of a CCDC segment, helps fill in gaps before the first segment, after the last segment, and between segments

# Tile/calval list options
clipping_method = 'tiles'  # Choose 'tiles' or 'calval'
skip_existing = True # Choose True = subset tile list to exclude tiles that have already been exported to the current CCDC image collection

# Reprojection options
scale = params.scale

# Output data options
in_version = 'v20240207'   # Input version identifier
out_version = 'v20240207'  # Output version identifier
out_dir = 'projects/arctic-biomass-mapping/assets/seasonal_percentile_doys/'
in_dir = 'projects/arctic-biomass-mapping/assets/CCDC_tiles/'
ccdc_path = in_dir + 'CCDC_C2_SR_' + clipping_method + '_' + str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + in_version
doy_path = in_dir + 'CCDC_C2_SR_' + clipping_method + '_' + str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + in_version + '_DOYs'
out_path = out_dir + 'seasonal_doys_' + clipping_method + '_' + str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + out_version

# 1.1 ----- EXISTING MODULES -----

from utils import misc
from utils import temporal_segmentation as temporal_seg
Segments = temporal_seg.Segments

# 1.2 ----- LOAD DATA -----

tiles_mgrs_s2 = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/tiles_mgrs_s2')
ccdc_fits = ee.ImageCollection(ccdc_path)
doy_fits = ee.ImageCollection(doy_path)

# ========================================
# 2. CREATE COLLECTIONS FOR PROCESSING ===
# ========================================

# 2.0 ----- CREATE COLLECTION OF TILES ACROSS WHICH TO MODEL REFLECTANCE -----

# Create image collection if it doesn't already exist
misc.createCollectionIfNotExists(out_path)

# If specified, remove tiles that have been already exported to the CCDC image collection
if skip_existing:
    existing_seasonal_doys = ee.ImageCollection(out_path).aggregate_array('system:index')
    ccdc_fits = ccdc_fits.filter(ee.Filter.inList('system:index', existing_seasonal_doys).Not())
    print('Existing tiles filtered out')

print('Tiles to export', ccdc_fits.aggregate_array('system:index').getInfo())
print('Final number of tiles:', ccdc_fits.aggregate_array('system:index').length().getInfo())

# 2.1 ----- CREATE LIST OF DATES ACROSS WHICH TO MODEL REFLECTANCE -----

# Create function to use input start/end year and start/end month/day to produce a list of dates
# Result is a list of lists:
#  - outer list = years
#  - inner lists = dates within years
def create_datelist(year):
    date_start = ee.Number(year).format('%.0f').cat('-').cat(start_MD)
    date_end = ee.Number(year).format('%.0f').cat('-').cat(end_MD)
    n_days = ee.Date(date_end).difference(ee.Date(date_start), 'day').round()
    days = ee.List.sequence(0, n_days, 1)

    def calc_dates(n):
        return ee.Date(date_start).advance(n, 'day')

    dates = days.map(calc_dates)
    return dates

date_list = map_years.map(create_datelist)
print('The list of years over which to model reflectance is:', map_years.getInfo())
print('The list of dates over which to model reflectance is:', date_list.getInfo())

# ==================================
# 3. CREATE PROCESSING FUNCTIONS ===
# ==================================

# 3.1 ----- SEASONAL DOY FUNCTION -----
# Create function to loop CCDC fit images (each corresponds to a MGRS tile) and:
# apply peak summer composite function to determine peak summer DOY for each year,
# composite annual peak summer DOYs to get final peak summer DOY,
# get start and end snow-free DOYs (calculated using Python script 3_0_proc_ccdc),
# and calculate early and late summer DOYs

# MAP OVER CCDC FIT IMAGES (EACH CORRESPONDS TO A MGRS TILE)
def season_doys_tiles(ccdc_img):
    
    # Get current tile name
    tile_name = ee.String(ccdc_img.id())

    # Get footprint from CCDC fit image
    footprint = ccdc_img.geometry()

    # Create CCDC segments for current CCDC fit image (corresponds to a MGRS tile)
    # Arranges the data to make it easier to work with
    # Attributes several functions to the resulting object
    ccdc_segments = Segments(ccdc_img, date_format=1)  # 1 = provide dates as fractional years

    # For each year's list of dates, find the max greenness and associated DOY.
    # ccdc_segments is captured from the enclosing scope via Python closure —
    # no wrapper function needed (unlike in JS GEE where closures aren't supported).
    def peak_summer_for_year(list_item):

        def map_dates(date):

            # Grab segment to use for peak summer calculation
            segment = ccdc_segments.find_by_date(ee.Date(date), strategy=segment_find_strategy)  # Expects ee.Date object

            # Produce modeled reflectance for current date
            fit = segment.slice(harmonics=3, extrapolateMaxDays=extrapolate_max_days, date=ee.Date(date))  # t: use fractional date, date: use ee.Date object

            # Convert date to day of year
            doy = ee.Date(date).getRelative('day', 'year').add(1)

            # Add day-of-year band
            # Mask 'doy' band to match fit object -- not sure if this actually reduces memory use...?
            fit = fit.addBands(ee.Image.constant(ee.Number(doy)).updateMask(fit.select(0).mask()).toUint16().rename('doy'))
            
            return fit
        
        fits = ee.ImageCollection(ee.List(list_item).map(map_dates))

        def select_index_doy(img):

            return img.select([peak_summer_index, 'doy'])
        
        return fits.map(select_index_doy).qualityMosaic(peak_summer_index)

    # Get maximum index and associated DOY across all years, one image per year
    peak_summer_ic = ee.ImageCollection(date_list.map(peak_summer_for_year))

    # Calculate median of annual max index composites
    peak_summer = peak_summer_ic.median().int16()

    # Select DOY only, rename to indicate seasonal percentile
    peak_summer = peak_summer.select('doy').rename('doy_peakSummer')

    # Get start and end of snowfree
    doy_img = doy_fits.filter(ee.Filter.eq('system:index', tile_name.replace('$', '_DOYs'))).first()
    start_snowfree = doy_img.select(start_snowfree_name).rename('doy_startSnowfree')
    end_snowfree = doy_img.select(end_snowfree_name).rename('doy_endSnowfree')

    # Break into two periods: start of snowfree season to peak, peak to end of snowfree season
    start_end_periods = peak_summer.subtract(start_snowfree).rename('days_startToPeak') \
        .addBands(end_snowfree.subtract(peak_summer).rename('days_peakToEnd'))

    # Calculate early summer (halfway between start of snowfree season and peak summer) and late summer (halfway between peak summer and end of snowfree season)
    doy_seasons = start_snowfree.add(start_end_periods.select('days_startToPeak').divide(2)).rename('earlySummer') \
        .addBands(peak_summer.add(start_end_periods.select('days_peakToEnd').divide(2)).rename('lateSummer'))
    early_summer = doy_seasons.select('earlySummer').rename('doy_earlySummer')
    late_summer = doy_seasons.select('lateSummer').rename('doy_lateSummer')

    # Combine
    seasonal_doys = start_snowfree.addBands(early_summer) \
                                  .addBands(peak_summer) \
                                  .addBands(late_summer) \
                                  .addBands(end_snowfree)

    return seasonal_doys.set('system:footprint', footprint) \
                        .set('id', tile_name.replace(in_version, out_version)) \
                        .uint16()

# =============================================
# 4. PROCESS TO PRODUCE SEASONAL DOY IMAGES ===
# =============================================

# Subset CCDC fit images to only index
# Confirmed that modeled fits are identical to fits produced across full CCDC fit images
def subset_ccdc_bands(img):
    img_metadata = img.select(['tStart', 'tEnd', 'tBreak', 'numObs', 'changeProb'])
    img_greenness = img.select('.*' + peak_summer_index + '.*')
    return img_metadata.addBands(img_greenness)

ccdc_fits = ccdc_fits.map(subset_ccdc_bands)

# Map over CCDC fit images
seasonal_doys_final = ccdc_fits.map(season_doys_tiles)

# =========================
# 5. EXPORT ===============
# =========================
# https://gis.stackexchange.com/questions/412094/computedobject-error-user-memory-limit-exceeded-when-exporting-small-extract

# Get list of image IDs
export_ids = seasonal_doys_final \
    .aggregate_array('system:index') \
    .getInfo()

# Pre-fetch tile name → EPSG mapping as a Python dict before the export loop.
# Each feature in tiles_mgrs_s2 has 'name' and 'epsg' properties; this fetches all of them
# in one GEE call so we can look up EPSG codes by tile name without a round-trip per tile.
tile_epsg = {
    feat['properties']['name']: feat['properties']['epsg']
    for feat in tiles_mgrs_s2.select(['name', 'epsg']).getInfo()['features']
}

# Export each image to image collection in assets
for img_id in export_ids:
    # Derive tile name from image ID using Python string methods (no GEE call needed)
    tile_name = img_id.replace('CCDC_C2_' + str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_', '').replace('_' + in_version, '')
    img = seasonal_doys_final.filter(ee.Filter.eq('system:index', img_id)).first()

    task = ee.batch.Export.image.toAsset(
        image=img,
        description=img_id,
        assetId=out_dir + 'seasonal_doys_' + clipping_method + '_' + str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + out_version + '/' + img_id,
        region=img.geometry(),
        crs='EPSG:' + tile_epsg[tile_name],
        scale=scale,
        maxPixels=1e13
    )
    task.start()

raise Exception('stop')
