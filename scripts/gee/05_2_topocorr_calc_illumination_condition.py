#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""

DESCRIPTION: Calculate illumination condition using Landsat scene DOY and sun angle/elevation information + DEM

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:
    - Gets solar information from the Landsat scene with DOY closest to seasonal DOY image, on a pixel-by-pixel basis
    - Landsat data is filtered and masked, but does not need to be pre-processed further because we are only using DOY and sun angle/elevation properties
    - For each MGRS tile Landsat scenes are filtered to exclude scenes outside the range given by the minimum DOY - maximum DOY across the tile
        - This information is calculated in a GEE script and provided in the 'tiles_mgrs_s2_grouped_start_doy_45n_doy_min_max_v20240207' featureCollection
    - Landsat scenes with WRS row >=35 are filtered out 
        - These scenes are at very high latitudes and have different overpass times from other data at these latitudes
        - Scenes with row >=35 are rare, most data at high latitudes have low WRS rows and similar overpass times
        - The difference in overpass time causes problems with the sun angle/elevation - DOY relationships

TO-DO:

"""


# In[13]:


# ===========================
# 1. SET UP =================
# ===========================

# 1.0 ----- PACKAGES -----

import ee 
import importlib

ee.Initialize(project='arctic-biomass-mapping')

lsat_proc = importlib.import_module('landsat_proc_mod')
lsat_xcal_ccdc = importlib.import_module('00_fun_landsat')


# In[14]:


# 1.1 ----- PARAMETERS -----

# Data input parameters
ccdc_start_year = 1984
ccdc_end_year = 2023
version = 'v20240207'

# Landsat filtering parameters
maxCloudCoverLand = 60; # Filtering out mostly cloudy scenes can speed processing a lot
start_year = 1984
end_year = 2023
wrs_row = 35

# Landsat masking parameters
waterOccurrence = 90; # Specify occurrence percentage threshold for classifying water
waterBuffer = 0; # Specify amount by which to buffer water (in meters)

# Specify season to process
seasons = ['startSnowfree', 'earlySummer', 'peakSummer', 'lateSummer', 'endSnowfree']


# In[15]:


# 1.2 ----- READ IN DATA -----

# ROIs
tiles = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/tiles_mgrs_s2_grouped_start_doy_45n_doy_min_max_v20240207')
arctic_poly = ee.FeatureCollection("projects/arctic-biomass-mapping/assets/ROIs/arctic_oroarctic_coast_buffer_laea_final")

# DOY data
doy_tiles = ee.ImageCollection('projects/arctic-biomass-mapping/assets/seasonal_percentile_doys/seasonal_doys_tiles_'+ str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + version)
doy_calval = ee.ImageCollection('projects/arctic-biomass-mapping/assets/seasonal_percentile_doys/seasonal_doys_calval_'+ str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + version)
doy = doy_tiles.merge(doy_calval).mosaic()

# DEM
# Copernicus DEM, 30 m resolution
dem = ee.ImageCollection("COPERNICUS/DEM/GLO30")             .mosaic()             .setDefaultProjection('EPSG:3857',None,30)             .select('DEM')             .rename('elevation')

# Landsat Collection 2 image collections
l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")


# In[16]:


# ======================
# 2. TIDY DATA =========
# ======================

# 2.2 ----- SET UP PROCESSING LIST -----

# Function to get grouped MGRS tile names as a list
def get_name(feat):
    
    return ee.Feature(feat).getString('name')

# Get grouped MGRS tile names as a list
tile_list = tiles.toList(50).map(get_name)

# Remove tiles with no reflectance data
tile_list = tile_list.filter(ee.Filter.inList('item', ['sdoy69', 'sdoy72', 'sdoy74', 'sdoy77', 'sdoy79']).Not())
print('Total number of tiles across ROI', ee.List(tile_list).length().getInfo())

# Convert tile name list to javascript list to enable exports
tile_list_js = tile_list.getInfo()

print('Tiles to export', tile_list_js)


# In[17]:


# 2.2 ----- PROCESS DEM -----

# Smooth the DEM using a low-pass kernel
boxcar = ee.Kernel.square(radius = 3, units = 'pixels', normalize = True)
dem_s = ee.Image(dem).convolve(boxcar)

# Create terrain layers
slp = ee.Terrain.slope(dem_s)
slp_rad = ee.Terrain.slope(dem_s).multiply(3.14159265359).divide(180)
asp_rad = ee.Terrain.aspect(dem_s).multiply(3.14159265359).divide(180)


# In[18]:


# ==============================================
# 3. CALCULATE ILLUMINATION CONDITION ==========
# ==============================================

for tile_name in tile_list_js: # LOOP TILES
    
    # 3.1 ----- FETCH AND TIDY MGRS TILE -----

    print('Processing tile: ', tile_name)
    
    # Get tile feature using tile name
    tile = tiles.filterMetadata('name','equals', tile_name).first()
    
    # Get intersection of full ROI (pan-Arctic polygon) and the current grouped MGRS tile
    roi = ee.Feature(tile).intersection(arctic_poly.geometry())
       
    ic_all = []
    for season in seasons: # LOOP SEASONS
        
        # 3.2 ----- GET SEASONAL MIN/MAX DOYS -----

        print('Current season', season)
    
        # Get seasonal start and end DOY
        # These values are the min and max from the seasonal DOY composite
        # Subtract/add 2 days in case closest Landsat scene is slightly outside this range
        start_doy = ee.Number(tile.get('doy_' + season + '_min')).subtract(2)
        end_doy = ee.Number(tile.get('doy_' + season + '_max')).add(2)

        # 3.3 ----- PRE-PROCESS LANDSAT DATA -----
        
        # No need to do any reflectance processing/cross-calibration etc.
        # All we need is the solar information properties and the doy band
        
        # Merge Landsat collections
        l5to8 = l5.merge(l7).merge(l8)

        # Filter Landsat collections by region, date, cloud cover, and WRS row
        l5to8f = l5to8.filterBounds(roi.geometry())                       .filter(ee.Filter.calendarRange(start_year, end_year, 'year'))                       .filter(ee.Filter.calendarRange(start_doy, end_doy))                       .filter(ee.Filter.lt('CLOUD_COVER_LAND', maxCloudCoverLand))                       .filter(ee.Filter.lt('WRS_ROW', wrs_row))

        # Tidy - rename bands        
        l5to8t = l5to8f.map(lsat_proc.l4to8_c2_rename_bands)

        # Apply masks
        l5to8m = l5to8t.map(lsat_proc.l4to8_c2_qa_maskClouds)                        .map(lsat_proc.l4to8_c2_qa_maskSnow)                        .map(lsat_xcal_ccdc.maskWater_jrc_occurrence(waterOccurrence, waterBuffer))                        .map(lsat_xcal_ccdc.maskOceanIce)
        
        # Add DOY band
        l5to8doy = l5to8m.map(lsat_xcal_ccdc.add_doy)

        # Select band(s)
        l5to8doy = l5to8doy.select('doy')
        
        # Get seasonal DOY data
        doy_season = doy.select('doy_' + season)
        
        # 3.4 ----- GET SOLAR INFORMATION FROM LANDSAT SCENES -----
        
        # Function to add solar information to Landsat scenes:
        # Get solar position data and convert to constant images
        # Calculate difference between Landsat scene DOY and seasonal DOY image
        def ls_solar(img):

            # Extract image metadata about solar position and convert to constant images
            SZ_rad = ee.Image.constant(ee.Number(90).subtract(ee.Number(img.get('SUN_ELEVATION')))).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000)).float()
            SA_rad = ee.Image.constant(ee.Number(img.get('SUN_AZIMUTH')).multiply(3.14159265359).divide(180)).clip(img.geometry().buffer(10000)).float()

            # Calculate difference between image DOY and seasonal DOY
            # Multiply by -1 so that less difference = higher number and we can use qualityMosaic
            doy_diff = img.select('doy').subtract(doy_season).abs().multiply(-1)

            # Return image with sun elevation, sun azimuth and DOY bands
            return img.select('doy').rename('doy_' + season).addBands(doy_diff.rename('doy_diff_' + season)).addBands(SZ_rad.rename('sz_rad')).addBands(SA_rad.rename('sa_rad'))

        # Map function across Landsat scenes
        l5to8_sun = l5to8doy.map(ls_solar)
       
        # Get solar information from Landsat scene with DOY closest to seasonal DOY image, on a pixel-by-pixel basis
        l5to8_final = l5to8_sun.qualityMosaic('doy_diff_' + season)

        # Get number of Landsat scenes availabile for quality mosaic
        l5to8_count = l5to8_sun.select('doy_' + season).count().rename('pixel_count_' + season)
        
        # 3.5 ----- CALCULATE ILLUMINATION CONDITION (IC) -----

        # Slope IC component
        cosZ = l5to8_final.select('sz_rad').cos()
        cosS = slp_rad.cos()
        slope_illumination = cosS.expression("cosZ * cosS",
                                              {'cosZ': cosZ,
                                               'cosS': cosS.select('slope')})

        # Aspect IC component
        sinZ = l5to8_final.select('sz_rad').sin()
        sinS = slp_rad.sin()
        cosAziDiff = (l5to8_final.select('sa_rad').subtract(asp_rad)).cos()
        aspect_illumination = sinZ.expression("sinZ * sinS * cosAziDiff",
                                               {'sinZ': sinZ,
                                                'sinS': sinS,
                                                'cosAziDiff': cosAziDiff})

        # Full IC
        ic = slope_illumination.add(aspect_illumination).rename('IC_' + season)
        
        # Tidy and add metadata
        ic = ic.multiply(10000).int16()                .addBands(l5to8_count.uint16())                .addBands(l5to8_final.select('doy_diff_' + season).multiply(-1).uint8())                .addBands(cosZ.rename('cosZ_' + season).multiply(10000).int16())                .addBands(cosS.rename('cosS_' + season).multiply(10000).int16())

        # Add to seasons list
        ic_all.append(ic)
        
    # 3.6  ----- TIDY FINAL IMAGE COLLECTION ----------

    # Convert to image collection
    # Convert to multiband image
    # Remove leading numbers
    # Scale by 10,000 for storage purposes
    ic_all = ee.ImageCollection.fromImages(ic_all)                                .toBands()                                .regexpRename('^[0-9]_', '')

    # ===================
    # 4. EXPORT =========
    # ===================
    
    task_export = ee.batch.Export.image.toAsset(
        image = ic_all.clip(roi.geometry()),
        description = 'illumination_condition_' + version + '_' + tile_name,
        assetId = 'projects/arctic-biomass-mapping/assets/topographic_correction/illumination_condition_'+ str(ccdc_start_year) + '_' + str(ccdc_end_year) + '_' + version + '/' + tile_name,
        region = roi.geometry(),
        scale = 30,
        crs = 'EPSG:3571',
        maxPixels = 1e13
    )

    task_export.start()

    print(tile_name, ': task submitted')
    


# In[ ]:




