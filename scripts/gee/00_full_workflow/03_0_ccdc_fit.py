"""

DESCRIPTION: Run CCDC across MGRS tiles

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:

TO-DO:

"""

# ===========================
# 1. SET UP =================
# ===========================

# 1.0 ----- PACKAGES -----

import ee
from utils import params


try:
    ee.Initialize(project=params.ee_project)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=params.ee_project)
    
from utils import landsat_proc_mod as lsat_proc
from utils import landsat as lsat_xcal_ccdc
from utils import ccdc as ccdc_run
from utils import misc

# 1.1 ----- PARAMETERS -----

# Tile/calval list options
skip_existing = True # Choose True = subset tile list to exclude tiles that have already been exported to the current CCDC image collection
utm_zone = ee.List.sequence(1, 60) # UTM zone as a list of numbers (e.g. [6, 12, 55]) or ee.List.sequence(1, 60), optional, but can use to further subset exported tiles; ee.List.sequence(1, 60) = do not subset by UTM zone -- NEW VERSION BASED ON ZONE NUMBER
calval_region = '' # Region name to subset exported cal/val data; 'USA', 'Finland', 'Sweden', 'Russia', 'Norway', 'Greenland', 'Canada'; '' = do not subset by region

# Choose 'tiles' or 'calval'
# Need to run both
# - tiles runs full pan-arctic image
# - calval gets smaller areas for calval sites outside pan-arctic boundaries
clipping_method = 'tiles' 

# Calval processing options
calval_buffer = ee.Number(900).divide(2) # 900m x 900m window to be safe
data_version = 'v20240215'

# Landsat processing options
start_year = 1984
end_year = 2023
maxCloudCoverLand = 60 # Filtering out mostly cloudy scenes can speed processing a lot
wrs_row = 35 # Filter out high WRS rows, they do not belong in the Arctic
waterOccurrence = 90 # Specify occurrence percentage threshold for classifying water
waterBuffer = 0 # Specify amount by which to buffer water (in meters)

# CCDC options
bands_ccdc = ['blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2', 'EVI2b', 'NBR', 'NDMI', 'NDVI', 'NDWI'] # Bands to run through CCDC; In order to apply cross-calibration direclty to indices, will need to run CCDC on indices rather than calculating later; remove quality bands and surface temperature (holes in surface temperature images cause issues)
breakpointBands = ['green', 'red', 'NIR', 'SWIR1', 'SWIR2'] # Bands to use for breakpoint detection in CCDC, changed to match global asset settings, blue likely subject to atmospheric noise
tmaskBands = ['green', 'SWIR2'] # The name or index of the bands to use for iterative TMask cloud detection
dateFormat = 1 # Decimal year, human understandable format
maxIterations = 10000 # Changed to match global asset, maybe runs a lot faster?
lam = 20 # Lambda: default of 20 makes sense when reflectance is scaled 0-10000 (0.002 if scaled 0-1)

# DOY options
snow_free_percentiles = [3, 97] # Percentiles to use for defining the snow-free season i.e. 'spring' and 'fall'
snow_free_percentile_names = ['p003', 'p097'] # Percentile names

# Output data options
version = 'v20240207' # Version identifier
asset_dir = 'projects/arctic-biomass-mapping/assets/CCDC_tiles/' # Asset path for export
ic_path_ccdc = asset_dir + 'CCDC_C2_SR_' + clipping_method + '_'+ str(start_year) + '_' + str(end_year) + '_' + version
ic_path_doy = asset_dir + 'CCDC_C2_SR_' + clipping_method + '_'+ str(start_year) + '_' + str(end_year) + '_' + version + '_DOYs'

# 1.2 ----- READ IN LANDSAT DATA -----

# Landsat Collection 2 image collections
l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")

# 1.3 ----- READ IN AND TIDY ROI DATA-----

# Biome ROI
arctic_poly = ee.FeatureCollection("projects/arctic-biomass-mapping/assets/ROIs/arctic_oroarctic_coast_buffer_laea_final")

# Clipping method specific tidying
if clipping_method == "tiles":
    
    # Read in tiles - Arctic specific
    tiles_mgrs_s2 = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/tiles_mgrs_s2_arctic_final')

    # Set ROI to clip to
    roi_filter_clip = arctic_poly

    # Set tile list
    tiles_mgrs_roi = tiles_mgrs_s2

elif clipping_method == "calval":
# When clipping_method = 'calval' we are only processing the CCDC for calval data that falls outside official Pan Arctic ROI

    # Read in tiles - full list
    tiles_mgrs_s2 = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/tiles_mgrs_s2')

    # Cal/val data
    calval_data = ee.FeatureCollection("projects/arctic-biomass-mapping/assets/field_data/arctic_tundra_biomass_synthesis_plots_" + data_version)
    
    # Filter to only calval data outside official Pan Arctic ROI
    calval_data = calval_data.filter(ee.Filter.intersects(".geo", arctic_poly.first().geometry()).Not());

    # Buffer cal/val data points
    calval_data_buffer = calval_data.map(lambda feat: feat.buffer(calval_buffer))

    # Dissolve buffered cal/val data points
    calval_data_union = calval_data_buffer.union()

    # Set ROI to clip to
    roi_filter_clip = calval_data_union

    # Set tile list
    tiles_mgrs_roi = tiles_mgrs_s2.filterBounds(roi_filter_clip)

else:
    raise ValueError("Unknown clipping method")

# ===============================
# 2. INITIATE TILE LIST =========
# ===============================

# Filter further by UTM zone
tiles_mgrs_roi = tiles_mgrs_roi.filter(ee.Filter.inList('zone', utm_zone))

# Get MGRS tile names as a list
tile_list = tiles_mgrs_roi.aggregate_array('name')
print('Total number of tiles across ROI:', ee.List(tile_list).length().getInfo())

# ==================================================
# 3. INITIATE/CHECK CCDC IMAGE COLLECTIONS =========
# ==================================================

# Create image collections if they don't already exist
misc.createCollectionIfNotExists(ic_path_ccdc)
misc.createCollectionIfNotExists(ic_path_doy)

# If specified, remove tiles that have been already exported to the CCDC image collection
if skip_existing:
    ccdc = ee.ImageCollection(ic_path_ccdc)
    ccdc_tiles = ccdc.aggregate_array('system:index').map(lambda s: ee.String(s).replace(f'CCDC_C2_{start_year}_{end_year}_', '').replace(f'_{version}', '').replace('_DOYs', ''))
    tile_list = tile_list.filter(ee.Filter.inList('item', ccdc_tiles).Not())
    print('Existing tiles filtered out')

# Convert tile name list to javascript list to enable exports
tile_list_js = tile_list.getInfo()

print('Final number of tiles:', ee.List(tile_list).length().getInfo())
print('Tiles to export', tile_list_js)

# ======================
# 4. RUN CCDC ==========
# ======================

# Report parameters
print('CCDC model fits will be applied on all Landsat data from ', start_year, ' to ', end_year)
print('Landsat scenes are filtered based on a max land cloud cover of: ', maxCloudCoverLand)
print('Landsat scenes are filtered based on a WRS row less than: ', wrs_row)
print('Landsat scenes are masked using a water occurence threshold of: ', waterOccurrence)
print('Landsat water mask is buffered by ', waterBuffer, ' meters')
print('The following bands are run through the CCDC: ', bands_ccdc)
print('Snow-free dates are established using the following percentiles: ', snow_free_percentile_names)
print('The following bands are used as CCDC breakpoint bands: ', breakpointBands)
print('The following bands are used as CCDC tMask bands: ', tmaskBands)
print('CCDC model fitting uses max iterations: ', maxIterations)
print('CCDC model fitting uses lambda: ', lam)
print('Results will be output to the following directory: ', asset_dir)
print('Results will be output to image collection with version identifier: ', version)

print('List of tiles to map over:', tile_list_js)

# Loop to run CCDC
# CCDC processing function automatically pre-processes Landsat 5, 7 and 8 data
# extractCcdcFromTile(tile_coll, roi, version, bands_ccdc, breakpointBands, tmaskBands, dateFormat, maxIterations, lambda, start_year, end_year, maxCloudCoverLand, asset_path)
for tile_name in tile_list_js:
    
    ccdc_run.extractCcdcFromTile(
        tile_name = tile_name,
        tile_coll = tiles_mgrs_roi,
        clipping_method = clipping_method,
        roi = roi_filter_clip,
        snow_free_percentiles = snow_free_percentiles,
        snow_free_percentile_names = snow_free_percentile_names,
        version = version,
        bands_ccdc = bands_ccdc,
        breakpointBands = breakpointBands,
        tmaskBands = tmaskBands,
        dateFormat = dateFormat,
        maxIterations = maxIterations,
        lam = lam,
        start_year = start_year,
        end_year = end_year,
        maxCloudCoverLand = maxCloudCoverLand,
        wrs_row = wrs_row,
        waterOccurrence = waterOccurrence,
        waterBuffer = waterBuffer,
        asset_dir = asset_dir
        
    )

raise Exception('stop')

