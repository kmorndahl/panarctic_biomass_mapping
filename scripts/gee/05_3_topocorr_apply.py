#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""

DESCRIPTION: Apply topographic correction to modeled seasonal reflectance imagery

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:
    - No need to intersect with arctic polygon for reduceRegion ROI because reflectance image is already clipped/masked to arctic polygon
    - DEM is processed and slope is calculated in a separate GEE script to reduce compute demands and eliminate error outs
    
TO-DO:

"""


# In[1]:


# ===========================
# 1. SET UP =================
# ===========================

# 1.0 ----- PACKAGES -----

import ee 
import importlib
import numpy

ee.Initialize(project='arctic-biomass-mapping')

lsat_proc = importlib.import_module('landsat_proc_mod')
lsat_xcal_ccdc = importlib.import_module('00_fun_landsat')


# In[2]:


# 1.1 ----- PARAMETERS -----

# Data input parameters
year = 2000
version = 'v20240207'
start_year = 1984
end_year = 2023

# Processing parameters
tile_scale = 1

# Topographic correction parameters
tile_pixels = 11100000 # Approximate number of pixels in a full MGRS tile
frac = 0.1 # Proportion of pixels to use to determine reflectance/illumination condition linear relationship

# Specify season to process
seasons = ['startSnowfree', 'earlySummer', 'peakSummer', 'lateSummer', 'endSnowfree']

# Band specification parameters
bands_topo_corr = ['spectral_blue_', 'spectral_green_', 'spectral_red_', 'spectral_NIR_', 'spectral_SWIR1_', 'spectral_SWIR2_', 'spectral_EVI2b_']; # No need to correct normalized spectral indices


# In[3]:


# 1.2 ----- READ IN DATA -----

# ROIs
tiles = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/tiles_mgrs_s2_arctic_final')
arctic_poly = ee.FeatureCollection("projects/arctic-biomass-mapping/assets/ROIs/arctic_oroarctic_coast_buffer_laea_final")

# DOY data
doy = ee.ImageCollection('projects/arctic-biomass-mapping/assets/seasonal_percentile_doys/seasonal_doys_tiles_' + str(start_year) + '_' + str(end_year) + '_' + version).mosaic()

# Modeled reflectance data
refl = ee.Image('projects/arctic-biomass-mapping/assets/seasonal_modeled_reflectance/seasonal_modeled_reflectance_' + str(start_year) + '_' + str(end_year) + '_' + version + '_' + str(year))

# Copernicus DEM, 30 m resolution
# Read in pre-processed DEM to save compute
slp = ee.Image('projects/arctic-biomass-mapping/assets/predictors/copernicus_dem_slope_clipped')

# Illumination condition
ic = ee.ImageCollection('projects/arctic-biomass-mapping/assets/topographic_correction/illumination_condition_' + str(start_year) + '_' + str(end_year) + '_' + version).mosaic()


# In[7]:


# ======================
# 2. TIDY DATA =========
# ======================

# 2.1 ----- SET UP PROCESSING LIST -----

print('Total number of tiles across ROI', tiles.size().getInfo())

# Function to get MGRS tile names as a list
def get_name(feat):
    
    return ee.Feature(feat).getString('name')

# Get MGRS tile names as a list
tile_list = tiles.toList(10000).map(get_name)

# Convert tile name list to javascript list to enable exports
tile_list_js = tile_list.getInfo()

# Report
print('Tiles to export', tile_list_js)
print('Number of tiles to export', len(tile_list_js))


# In[8]:


# ====================================
# 3. TOPOGRAPHIC CORRECTION ==========
# ====================================

for tile_name in tile_list_js: # LOOP TILSE

    print('Processing tile: ', tile_name)

    # 3.1 ----- FETCH AND TIDY MGRS TILE -----

    # Get tile feature using tile name
    tile = tiles.filterMetadata('name','equals', tile_name).first()
          
    # Get pixel count for determining sample size for calculating reflectance/illumination condition linear relationship
    pixel_count = ee.Number(tile_pixels) # From individual MGRS tiles
   
    refl_topo_corr = []
    for season in seasons: #  LOOP SEASONS
        
        print('Current season', season)
        
        # ----- 3.2 PROCESS REFLECTANCE DATA -----
        
        # NOTE: confirmed calculations are valid using reflectance data scaled by 10,000

        # Get seasonal illumination condition (IC)
        ic_season = ic.select('IC_' + season).divide(10000)
        cosZ_season = ic.select('cosZ_' + season).divide(10000)
        cosS_season = ic.select('cosS_' + season).divide(10000)

        # Add IC
        refl_ic = ee.Image(refl.divide(10000)
                               .addBands(ic_season) \
                               .addBands(cosZ_season.rename('cosZ')) \
                               .addBands(cosS_season.rename('cosS')) \
                               .addBands(slp.rename('slope')))

        # Mask slope - no need to topographically correct flat areas
        slope_mask = refl_ic.select('slope').gte(5)
        refl_ic = ee.Image(refl_ic.updateMask(slope_mask))
        
        # Mask out IC values <= 0
        # These are very difficult to correct: https:#www.ingentaconnect.com/content/asprs/pers/2012/00000078/00000009/art00005?crawler=True
        ic_mask = refl_ic.select('IC_' + season).gt(0)
        refl_ic = ee.Image(refl_ic.updateMask(ic_mask))
        
        # Get bands for current season
        season_bands = []
        for band in bands_topo_corr:
            season_bands.append(band + season)
        print('Bands to topographically correct:')
        print(season_bands)
               
        # Count pixels - cannot topographically correct if too few unmasked pixels
        num_pixels = refl_ic.select('spectral_blue_' + season).reduceRegion(           reducer = ee.Reducer.count(),           geometry = tile.geometry(),           scale = 30,           maxPixels = 1e13,           tileScale = tile_scale) # Adjust tileScale as necessary to account for out of memory errors and computation timed out errors
        num_pixels = ee.Number(num_pixels.values().get(0))
               
        # ----- 3.3 APPLY TOPOGRAHIC CORRECTION -----
       
        refl_SCSccorr = []
        for band in season_bands: # LOOP BANDS
            
            # Linear fit using random sample of pixels
            random_pixels = refl_ic.select(ee.List(['IC_' + season, band]))                                    .addBands(refl_ic.select(band).mask().rename('class').byte())                                    .stratifiedSample(numPoints = pixel_count.multiply(frac).int(),                                                      classBand = 'class',                                                      region = tile.geometry(),                                                      scale = 30,                                                      seed = 42,                                                      dropNulls = True,                                                      geometries = False,                                                      tileScale = tile_scale); # Adjust tileScale as necessary to account for out of memory errors and computation timed out errors
            out = random_pixels.reduceColumns(reducer = ee.Reducer.linearFit(), selectors = ['IC_' + season, band])
                
            # Compute linear fit coefficients: a(slope), b(intercept), c(b/a)
            out_a = ee.Number(out.get('scale'))
            out_b = ee.Number(out.get('offset'))
            out_c = ee.Number(out.get('offset')).divide(ee.Number(out.get('scale')))

            # Apply the SCSc correction
            SCSc_output = refl_ic.expression("((image * (cosB * cosZ + cvalue)) / (ic + cvalue))", {
              'image': refl_ic.select(ee.List([band])),
              'ic': refl_ic.select('IC_' + season),
              'cosB': refl_ic.select('cosS'),
              'cosZ': refl_ic.select('cosZ'),
              'cvalue': out_c
              })

            # Fill masked areas with original reflectance values
            SCSc_output = SCSc_output.multiply(10000)                                      .unmask(refl.select(band), False); # Must choose False otherwise output is empty

            # Add to list of topographically corrected bands
            refl_SCSccorr.append(ee.Image(SCSc_output))
        
        # Convert to multiband image and remove leading numbers
        refl_SCSccorr = ee.ImageCollection.fromImages(refl_SCSccorr).toBands().regexpRename('^[0-9]_', '');

        # Add topographically corrected seasonal image to list of topographically corrected images
        # If fewer than 100 unmasked pixels available for determining illumination condition/reflectance relationship...
        # Do not perform topographic correction and instead return original (uncorrected) reflectance image
        refl_append = ee.Algorithms.If(num_pixels.lt(100), refl.select(season_bands), refl_SCSccorr)
        refl_topo_corr.append(refl_append)

    # 3.4  ----- TIDY FINAL IMAGE COLLECTION ----------

    # Convert to multiband image and remove leading numbers
    refl_topo_corr = ee.ImageCollection.fromImages(refl_topo_corr).toBands().regexpRename('^[0-9]_', '');

    # ===================
    # 4. EXPORT =========
    # ===================

    task_export = ee.batch.Export.image.toAsset(
        image = refl_topo_corr.int16(), # Convert to int 16 to save storage space
        description = 'seasonal_modeled_reflectance_' + str(start_year) + '_' + str(end_year) + '_' + version + '_' + str(year) + '_topographic_correction_' + tile_name,
        assetId = 'projects/arctic-biomass-mapping/assets/seasonal_modeled_reflectance/seasonal_modeled_reflectance_' + str(start_year) + '_' + str(end_year) + '_' + version + '_' + str(year) + '_topographic_correction/' + tile_name,
        region = tile.geometry(),
        scale = 30,
        crs = 'EPSG:3571',
        maxPixels = 1e13
    )

    task_export.start()

    print(tile_name, ': task submitted')
    


# In[ ]:




