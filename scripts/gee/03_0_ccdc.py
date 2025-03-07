#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""

DESCRIPTION: Run CCDC across MGRS tiles

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:

TO-DO:

"""


# In[7]:


# Only need to run once, credentials are stored
# https://developers.google.com/earth-engine/guides/python_install#authentication
# https://developers.google.com/earth-engine/guides/auth
ee.Authenticate()


# In[1]:


# ===========================
# 1. SET UP =================
# ===========================

# 1.0 ----- PACKAGES -----

import ee 
import importlib

ee.Initialize(project='arctic-biomass-mapping')

lsat_proc = importlib.import_module('landsat_proc_mod')
lsat_xcal_ccdc = importlib.import_module('00_fun_landsat')
ccdc_run = importlib.import_module('00_fun_ccdc')


# In[2]:


# 1.1 ----- PARAMETERS -----

# Tile/calval list options
append_ccdc = 'no' # Choose 'yes' = appending tiles/updating tiles for existing CCDC image collection; 'no' = creating new CCDC image collection
check_existing = True # Choose True = subset tile list to exclude tiles that have already been exported to the current CCDC image collection
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
wrs_row = 35
waterOccurrence = 90 # Specify occurrence percentage threshold for classifying water
waterBuffer = 0 # Specify amount by which to buffer water (in meters)
bands_tc = ['blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2', 'EVI2b']

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
asset_path = 'projects/arctic-biomass-mapping/assets/CCDC_tiles/' # Asset path for export


# In[3]:


# 1.2 ----- READ IN DATA -----

# Landsat Collection 2 image collections
l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")

# ROI
arctic_poly = ee.FeatureCollection("projects/arctic-biomass-mapping/assets/ROIs/arctic_oroarctic_coast_buffer_laea_final")

# Cal/val data
calval_data = ee.FeatureCollection("projects/arctic-biomass-mapping/assets/field_data/arctic_tundra_biomass_synthesis_plots_" + data_version) # Full dataset

# MGRS tiles
tiles_mgrs_s2 = ee.FeatureCollection(ee.Algorithms.If(ee.String(clipping_method).equals('calval'),                                                       ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/tiles_mgrs_s2'),                                                       ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/tiles_mgrs_s2_arctic_final')));
new_tiles = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/updated_dataset_append_ccdc/new_ds_tiles_20230601') # If adding to existing CCDC output, specify new tiles

# Copernicus DEM, 30 m resolution
dem = ee.ImageCollection("COPERNICUS/DEM/GLO30")             .mosaic()             .setDefaultProjection('EPSG:3857',None,30)             .select('DEM')             .rename('elevation')


# In[4]:


# ======================
# 2. TIDY DATA =========
# ======================

# 2.1 ----- PROCESS DEM -----

# Smooth the DEM using a low-pass kernel
boxcar = ee.Kernel.square(radius = 3,
                          units = 'pixels',
                          normalize = True)
dem_s = ee.Image(dem).convolve(boxcar)


# In[5]:


# 2.2 ----- TIDY CAL/VAL DATA -----
# When clipping_method = 'calval' we are only processing the CCDC for calval data that falls outside official Pan Arctic ROI

# Filter to only calval data outside official Pan Arctic ROI
calval_data = calval_data.filter(ee.Filter.intersects(".geo", arctic_poly.first().geometry()).Not());

# Function to buffer cal/val data points
def buffer_pts(feat):  
    
    return feat.buffer(calval_buffer)

# Buffer cal/val data points
calval_data_buffer = calval_data.map(buffer_pts)

# Dissolve buffered cal/val data points
calval_data_union = calval_data_buffer.union()

# Set ROI
roi_filter_clip = ee.FeatureCollection(ee.Algorithms.If(ee.String(clipping_method).equals('calval'),                                                         calval_data_union,                                                         arctic_poly))


# In[6]:


# 2.3 ----- SET UP PROCESSING LIST -----

# Get MGRS tiles for study area
tiles_mgrs_roi = ee.FeatureCollection(ee.Algorithms.If(ee.String(clipping_method).equals('calval'),                                                        tiles_mgrs_s2.filterBounds(roi_filter_clip),                                                        tiles_mgrs_s2));

# Filter further by UTM zone if specified
tiles_mgrs_roi = tiles_mgrs_roi.filter(ee.Filter.inList('zone', utm_zone))

# Function to get MGRS tile names as a list
def get_name(feat):
    
    return ee.Feature(feat).getString('name')

# Get MGRS tile names as a list
tile_list = tiles_mgrs_roi.toList(5000).map(get_name)
print('Total number of tiles across ROI', ee.List(tile_list).length().getInfo())

# Function to rename CCDC and DOY collection IDs to match tiles
def fix_tile_name(string):
    return ee.String(string).replace('CCDC_C2_' + str(start_year) + '_' + str(end_year) + '_', '').replace('_' + version, '').replace('_DOYs', '')

# If specified, remove tiles that have been already exported to the CCDC image collection
if check_existing:
    doy = ee.ImageCollection('projects/arctic-biomass-mapping/assets/CCDC_tiles/CCDC_C2_SR_' + clipping_method + '_' + str(start_year) + '_' + str(end_year) + '_' + version + '_DOYs')
    ccdc = ee.ImageCollection('projects/arctic-biomass-mapping/assets/CCDC_tiles/CCDC_C2_SR_' + clipping_method + '_' + str(start_year) + '_' + str(end_year) + '_' + version)
    doy_tiles = doy.aggregate_array('system:index').map(fix_tile_name)
    ccdc_tiles = ccdc.aggregate_array('system:index').map(fix_tile_name)
    tile_list = tile_list.filter(ee.Filter.inList('item', ccdc_tiles).Not())
    print('Existing CCDC tile list size', ee.List(ccdc_tiles).length().getInfo())
    print('Existing DOY tile list size', ee.List(doy_tiles).length().getInfo())

# Finalize tile list
tile_list = ee.Algorithms.If(append_ccdc == 'yes', new_tiles.aggregate_array('name'), tile_list)
print('Final tile list size', ee.List(tile_list).length().getInfo())

# Convert tile name list to javascript list to enable exports
tile_list_js = tile_list.getInfo()

print('Tiles to export', tile_list_js)


# In[7]:


# ======================
# 3. RUN CCDC ==========
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
print('Results will be output to the following folder: ', asset_path)
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
        dem = dem_s,
        bands_ccdc = bands_ccdc,
        breakpointBands = breakpointBands,
        tmaskBands = tmaskBands,
        bands_tc = bands_tc,
        dateFormat = dateFormat,
        maxIterations = maxIterations,
        lam = lam,
        start_year = start_year,
        end_year = end_year,
        maxCloudCoverLand = maxCloudCoverLand,
        wrs_row = wrs_row,
        waterOccurrence = waterOccurrence,
        waterBuffer = waterBuffer,
        asset_path = asset_path
        
    )


# In[ ]:




