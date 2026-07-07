"""

DESCRIPTION: Run Monte Carlo regression forests (fit in R) across Pan Arctic domain

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:
- Reprojection of individual thematic predictor sets NOT necessary
- Predictors tried but removed from modeling:
  - MODIS year since fire:
    - The database only goes back to 2001... so for areas that haven't burned in this period, what fire year to assign? Could have burned in 2000 so fire year assignment might be misleading
    - Some calval data has unknown harvest day so if there was a fire in the same year it would be hard to determine if it occurred before or after harvest
      - This is probably the least concerning issue because presumably folks aren't sampling biomass in an area that burned earlier that same year (or if so it would be noted)
    - Same issue for generating wall-to-wall predictors for 2000 and 2020, if a fire occured in 2000 or 2020 difficult to know what fire year to assign (0? or next most recent fire? depends on when in year fire occured)
    - NOTE: new burned area product goes back to 1985 (but is missing some years...): https://gee-community-catalog.org/projects/gabam/#data-preprocessing
  - ESA permafrost data
    - Lots of sporadic data gaps
    - Large data gap in NE Greenland
- Exporting global images without gaps:
  - https://groups.google.com/g/google-earth-engine-developers/c/55btJBTGQoI/m/qXk0ce0XEgAJ
  - https://groups.google.com/g/google-earth-engine-developers/c/cUm_4WkiiAQ/m/6tQcuxLcBgAJ

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

# 1.0 ----- PARAMETERS -----

# CCDC date range
start_year = 1984
end_year = 2023

# Model options
ds_type = 'woody'
year = 2000  # Select year for predictions
date = ee.Date.fromYMD(year, 7, 31)  # Create date for extracting segment, choose mid-summer-ish for month and day

# Reprojection options
crs = params.crs
scale = params.scale  # Scale

# Output data options
ccdc_version = 'v20240207'
model_version = 'v20240514'
clipping_method = 'tiles'  # Choose 'tiles' or 'calval'
asset_path = 'projects/arctic-biomass-mapping/assets/modeled_continuous/model_continuous_' + str(year) + '_' + ds_type + '_' + model_version + '_mc/'  # Asset path for export
mc_list = ee.List.sequence(91, 100, 1).getInfo()

print('Monte carlo iterations: ', mc_list)

# 1.1 ----- EXISTING MODULES/FUNCTIONS -----

from utils import misc as fun_misc
from utils import temporal_segmentation as temporal_seg
Segments = temporal_seg.Segments

# Ensure output ImageCollection exists
fun_misc.createCollectionIfNotExists(asset_path.rstrip('/'))

# ======================================
# 2. MAP OVER MONTE CARLO ITERATIONS ===
# ======================================

for mc in mc_list:

    print('===== ITERATION  =====', mc)

    # 2.0 ----- RETRIEVE MODEL -----

    # Get model string from Google Storage
    mod_str = ee.List([ee.Blob('gs://arctic_biomass_mapping_models/' + model_version + '/mc/continuous/final_trees/' + ds_type + '_continuous_formatted_gee_short_' + model_version + '_' + str(mc) + '.txt').string()])
    print('Model string:', mod_str)
    print('Model string length: ', mod_str.join("\n").length())

    # Convert to GEE classifier
    mod = ee.Classifier.decisionTreeEnsemble(mod_str).setOutputMode('REGRESSION')
    print('Model:', mod)

    # Get information about the trained classifier
    print('Results of trained classifier', mod.explain())

    # 2.2 ----- EXISTING DATA PATHS -----

    # CCDC
    phenology = ee.ImageCollection('projects/arctic-biomass-mapping/assets/seasonal_percentile_doys/seasonal_doys_tiles_' + str(start_year) + '_' + str(end_year) + '_' + ccdc_version).mosaic()
    refl = ee.Image('projects/arctic-biomass-mapping/assets/seasonal_modeled_reflectance/seasonal_modeled_reflectance_' + str(start_year) + '_' + str(end_year) + '_' + ccdc_version + '_' + str(year))
    refl_tc = ee.ImageCollection('projects/arctic-biomass-mapping/assets/seasonal_modeled_reflectance/seasonal_modeled_reflectance_' + str(start_year) + '_' + str(end_year) + '_' + ccdc_version + '_' + str(year) + '_topographic_correction')

    # Model information
    selected_predictors = ee.FeatureCollection("projects/arctic-biomass-mapping/assets/model_ready_data/" + model_version + "/mc/continuous/final_predictors/final_predictors_" + ds_type + "_continuous_" + model_version + '_' + str(mc))

    # 2.3 ----- TIDY DATA -----

    # Reflectance data
    refl_tc = refl_tc.mosaic()  # Mosaic topographically corrected collection
    refl = refl.select('.*_ND.*|.*NBR.*')  # Get normalized indices from original image
    modeled_reflectance = refl_tc.addBands(refl)  # Combine

    # Model data
    def replace_dots(item):
        return ee.String(item).replace('[.]', '', 'g')

    def get_category(predictor_name):
        return ee.String(predictor_name).replace('_.*$', '')

    selected_predictors = ee.List(selected_predictors.aggregate_array('predictor').map(replace_dots))
    predictor_categories = selected_predictors.map(get_category)
    print('Selected predictors: ', selected_predictors)
    print('Predictor categories: ', predictor_categories)

    # ================================================================
    # 3. PREDICTORS ==================================================
    # ================================================================

    # 3.0 ----- LATITUDE/LONGITUDE PREDICTORS -----

    lat_lon = ee.Image.pixelLonLat()

    #  3.1 ----- TEXTURE PREDICTORS -----

    texture = ee.Image("projects/arctic-biomass-mapping/assets/predictors/seasonal_modeled_reflectance_NDVI_texture_" + str(year) + "_" + ccdc_version)
    print('Texture:', texture)

    #  3.2 ----- TOPOGRAPHIC PREDICTORS -----

    # MERIT DEM based, 90m resolution -----
    # https://gee-community-catalog.org/projects/geomorpho90/

    # Import
    dem = ee.Image("MERIT/DEM/v1_0_3").rename('topo_dem')
    cti = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/cti").mosaic().rename('topo_cti')
    slope = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/slope").mosaic().rename('topo_slope')
    tpi = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/tpi").mosaic().rename('topo_tpi')
    eastness = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/eastness").mosaic().rename('topo_eastness')
    northness = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/northness").mosaic().rename('topo_northness')

    # Topographic predictors
    topographic = (dem.addBands(cti)
                      .addBands(slope)
                      .addBands(tpi)
                      .addBands(eastness)
                      .addBands(northness))

    topographic_fill_gaps = topographic.focalMean(radius=dem.projection().nominalScale().multiply(3).divide(2), kernelType='square', units='meters', iterations=1)  # Fill some gaps to avoid NAs in cal/val data
    topographic = topographic_fill_gaps.blend(topographic)  # Fill some gaps to avoid NAs in cal/val data

    print('Topographic predictors:', topographic)

    #  3.3 ----- PERMAFROST PREDICTORS -----

    # Gruber permafrost data
    # https://tc.copernicus.org/articles/6/221/2012/tc-6-221-2012.pdf
    permafrost = ee.Image('projects/arctic-biomass-mapping/assets/predictors/global_permafrost_zonation_index_gruber2012_filled_20km_v3').rename('permafrost_index')

    print('Permafrost:', permafrost)

    #  3.4 ----- LANDCOVER PREDICTORS -----

    # Get World Terrestrial Ecosystems data
    # https://www.sciencebase.gov/catalog/item/6296791ed34ec53d276bb293
    # https://rmgsc.cr.usgs.gov/outgoing/ecosystems/Global/
    wte = ee.Image('projects/arctic-biomass-mapping/assets/predictors/wte_2020_filled_encoded')
    wte = wte.updateMask(wte.mask().gt(0))  # Remove mask percentage

    print('World Terrestrial Ecosystems:', wte)

    # Ecoregions -----

    ecoregions = ee.Image('projects/arctic-biomass-mapping/assets/predictors/ecoregions_img')

    print('Ecoregions:', ecoregions)

    # Arctic zones -----

    zones_final = ee.Image('projects/arctic-biomass-mapping/assets/predictors/zones_img')

    print('Zones:', zones_final)

    #  3.5 ----- TREE COVER PREDICTORS -----

    # Get Hansen tree cover data
    tree_cover = ee.Algorithms.If(
        ee.Number(year).eq(2020),
        ee.Image("projects/arctic-biomass-mapping/assets/predictors/tree_cover_2020"),
        (ee.Image("UMD/hansen/global_forest_change_2022_v1_10")
           .select(['treecover2000'], ['cover'])  # Select relevant bands and rename
           .regexpRename('^', 'trees_')  # Add prefix
           .addBands(ee.Image("UMD/hansen/global_forest_change_2022_v1_10").select('treecover2000').gt(0).rename('trees_presence'))  # Calculate presence/absence
           .unmask(0))  # Unmask to set unmapped high latitude and ocean areas to 0 cover and 0 presence
    )

    print('Tree cover:', tree_cover)

    # 3.6 ----- MODELED REFLECTANCE AUXILIARY PREDICTORS -----

    # Calculate annual spectral predictors and rates of change
    # Per-band summaries (mean, median, range); change rates batched separately below
    # to keep the GEE graph compact (see calcAllChangeRates)
    spectral = (fun_misc.calcAnnualSpectral(modeled_reflectance, 'blue')
                .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'green'))
                .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'red'))
                .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'NIR'))
                .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'SWIR1'))
                .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'SWIR2'))
                .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'EVI2b'))
                .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'NBR'))
                .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'NDMI'))
                .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'NDVI'))
                .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'NDWI'))
                .addBands(fun_misc.calcAllChangeRates(modeled_reflectance, phenology)))

    print('Spectral:', spectral)

    # 3.7 ----- COMBINE PREDICTORS -----

    # Combine predictors
    predictors = (spectral.addBands(lat_lon)
                          .addBands(permafrost)
                          .addBands(texture)
                          .addBands(wte.reduceResolution(reducer=ee.Reducer.mode(), maxPixels=1000))
                          .addBands(ee.Image.constant(0).rename('world_terrestrial_ecosystems_other').setDefaultProjection(crs=crs, scale=scale).reduceResolution(reducer=ee.Reducer.mode(), maxPixels=1000))
                          .addBands(tree_cover)
                          .addBands(zones_final.reduceResolution(reducer=ee.Reducer.mode(), maxPixels=1000))
                          .addBands(ecoregions.reduceResolution(reducer=ee.Reducer.mode(), maxPixels=1000)))

    # Topographic predictors require computationally expensive focal mean, only add if necessary
    predictors = ee.Image(ee.Algorithms.If(predictor_categories.contains('topo'),
                                           predictors.addBands(topographic),
                                           predictors))

    print('Final predictors: ', predictors)

    # NOTE: checked that none of the predictors have a gap in E Russia

    # ===========================================================
    # 4. MODEL ==================================================
    # ===========================================================

    # Apply model
    results = predictors.classify(mod, 'predicted')

    # Return to appropriate scale
    results = results.round().uint16()
    print('Final results:', results)

    # ============================================================
    # 5. EXPORT ==================================================
    # ============================================================

    task = ee.batch.Export.image.toAsset(
        image=results,
        description='model_continuous_' + str(year) + '_' + ds_type + '_' + model_version + '_' + str(mc),
        assetId=asset_path + str(mc),
        region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
        scale=scale,
        crs=crs,
        maxPixels=1e13
    )
    task.start()

raise Exception('stop')