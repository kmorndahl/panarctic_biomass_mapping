import ee
from . import landsat_proc_mod as lsat_proc
from . import landsat as lsat_xcal_ccdc

# FUNCTION: format_ls_imagery
# USE: Scale Landsat reflectance to store in 16-bit integer, copy properties
# PARAMETERS:
# NOTES:
# AUTHOR: Kathleen Orndahl
# LAST UPDATE: 11-10-2024
# TO-DO:

def format_ls_imagery(img):

    result = img.multiply(10000).int16()
    result = result.copyProperties(img) \
                   .copyProperties(img, ['system:footprint']) \
                   .copyProperties(img, ['system:index']) \
                   .copyProperties(img, ['system:time_end']) \
                   .copyProperties(img, ['system:time_start'])

    return result

# FUNCTION: extractCcdcFromTile
# USE: Create CCDC model fits from input Landsat data for a regional tile
# PARAMETERS:
# NOTES:
#  - Loop over list of MGRS tile names
#  - Uses function wrapping to allow multiple parameters in .map()
#  - https://gis.stackexchange.com/questions/302760/gee-imagecollection-map-with-multiple-input-function
#  - https://stackoverflow.com/questions/56195985/can-i-map-a-function-with-more-then-one-argument
# AUTHOR: Matt Macander
# LAST UPDATE:
# TO-DO:

def extractCcdcFromTile(**kwargs):

    # Specify named arguments
    tile_name = kwargs['tile_name']
    tile_coll = kwargs['tile_coll'] # Feature collection to use for tiling
    clipping_method = kwargs['clipping_method'] # Method for clipping output, either by MGRS tiles or by buffered cal/val plots
    roi = kwargs['roi'] # Region of interest, used for clipping
    snow_free_percentiles = kwargs['snow_free_percentiles'] # Percentiles to use for defining the snow-free season i.e. 'spring' and 'fall'
    snow_free_percentile_names = kwargs['snow_free_percentile_names'] # Percentile names
    version = kwargs['version'] # Version identifier
    bands_ccdc = kwargs['bands_ccdc'] # Bands to run through CCDC'] In order to apply cross-calibration direclty to indices, will need to                                           run CCDC on indices rather than calculating later'] remove quality bands and surface temperature                                             (holes in surface temperature images cause issues)
    breakpointBands = kwargs['breakpointBands'] # Bands to use for breakpoint detection in CCDC
    tmaskBands = kwargs['tmaskBands'] # The name or index of the bands to use for iterative TMask cloud detection. These are typically the                                           green band and the SWIR1 band.
    dateFormat = kwargs['dateFormat'] # The time representation to use during fitting': 0 = jDays, 1 = fractional years, 2 = unix time in                                           milliseconds.
    maxIterations = kwargs['maxIterations'] # Maximum number of runs for LASSO regression convergence
    lam = kwargs['lam'] # Lambda for LASSO regression fitting
    start_year = kwargs['start_year'] # Start year
    end_year = kwargs['end_year'] # End year
    maxCloudCoverLand = kwargs['maxCloudCoverLand'] # Maximum cloud cover used during Landsat pre-processing
    wrs_row = kwargs['wrs_row'] # Maximum WRS row used during Landsat pre-processing
    waterOccurrence = kwargs['waterOccurrence'] # Threshold for masking out water, occurrence percentage = what percentage of mapped years                                                     water is presesnt
    waterBuffer = kwargs['waterBuffer'] # Buffer to apply to water mask, in meters, 0 = no buffer
    asset_dir = kwargs['asset_dir'] # Full path for where to export asset, include trailing slash

    # PROCESS MGRS TILE INFORMATION -----

    # Get tile feature using tile name
    footprint = tile_coll.filterMetadata('name','equals', tile_name).first()

    # Clip to exclude ocean, but include full tile area for tiles along edge of larger study area
    # For calval data, clip to buffered plot/site points
    clipper = ee.Feature(footprint).intersection(roi.geometry())

    # Get EPSG
    epsg = footprint.getString('epsg').getInfo()

    # Start and end DOY based on latitude, stored in tile columns
    start_doy = footprint.getNumber('start_doy')
    end_doy = footprint.getNumber('end_doy')

    # LANDSAT PRE-PROCESSING -----

    # Load Landsat Collection 2 SR image collections
    l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
    l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
    l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")

    # Filter Landsat SR collections by region, year, tile specific day-of-year, cloud cover and WRS row
    l5f = l5.filterBounds(clipper.geometry()) \
            .filter(ee.Filter.calendarRange(start_doy, end_doy)) \
            .filter(ee.Filter.calendarRange(start_year, end_year, 'year')) \
            .filter(ee.Filter.lt('CLOUD_COVER_LAND', maxCloudCoverLand)) \
            .filter(ee.Filter.lt('WRS_ROW', wrs_row))
    l7f = l7.filterBounds(clipper.geometry()) \
            .filter(ee.Filter.calendarRange(start_doy, end_doy)) \
            .filter(ee.Filter.calendarRange(start_year, end_year, 'year')) \
            .filter(ee.Filter.lt('CLOUD_COVER_LAND', maxCloudCoverLand)) \
            .filter(ee.Filter.lt('WRS_ROW', wrs_row))
    l8f = l8.filterBounds(clipper.geometry()) \
            .filter(ee.Filter.calendarRange(start_doy, end_doy)) \
            .filter(ee.Filter.calendarRange(start_year, end_year, 'year')) \
            .filter(ee.Filter.lt('CLOUD_COVER_LAND', maxCloudCoverLand)) \
            .filter(ee.Filter.lt('WRS_ROW', wrs_row))

    # Tidy Landsat collections for additional pre-processing:
    # Rename bands
    # Apply scale and offset
    # Calculate indices
    # Select bands
    l5t = l5f.map(lsat_proc.l4to8_c2_rename_bands) \
             .map(lsat_proc.l4to8_c2_scaleoff) \
             .map(lsat_xcal_ccdc.l4to8_c2_indices_xcal)
    l7t = l7f.map(lsat_proc.l4to8_c2_rename_bands) \
             .map(lsat_proc.l4to8_c2_scaleoff) \
             .map(lsat_xcal_ccdc.l4to8_c2_indices_xcal)
    l8t = l8f.map(lsat_proc.l4to8_c2_rename_bands) \
             .map(lsat_proc.l4to8_c2_scaleoff) \
             .map(lsat_xcal_ccdc.l4to8_c2_indices_xcal)

    # Cross-calibrate Landsat 5 and 8 to match Landsat 7
    # Use Arctic specific calibration parameters
    l5x = l5t.map(lsat_xcal_ccdc.applyCrossCalibrationLS5)
    l7x = l7t
    l8x = l8t.map(lsat_xcal_ccdc.applyCrossCalibrationLS8)

    # Merge collections
    l5to8 = l5x.merge(l7x).merge(l8x)

    # Mask clouds
    # Mask snow
    # Mask water
    # Mask ocean and Greenland icesheet
    # Add day of year band
    l5to8m = l5to8.map(lsat_proc.l4to8_c2_qa_maskClouds) \
                  .map(lsat_proc.l4to8_c2_qa_maskSnow) \
                  .map(lsat_xcal_ccdc.maskWater_jrc_occurrence(waterOccurrence, waterBuffer)) \
                  .map(lsat_xcal_ccdc.maskOceanIce) \
                  .map(lsat_xcal_ccdc.add_doy)

    # Check for valid pixels
    if(l5to8m.size().getInfo() == 0):
        print(tile_name, ': no valid pixels found, skipping tile export')
        return None

    # CREATE DAY-OF-YEAR PERCENTILE COMPOSITES -----

    # Select the day of year band and apply percentile reducer
    # For each pixel gets the day of year associated with each of the percentiles
    doys = l5to8m.select(['doy']) \
                 .reduce(ee.Reducer.percentile(snow_free_percentiles, snow_free_percentile_names)) \
                 .set('version', version)

    # RUN CCDC -----

    # Select the bands to use for the CCDC algorithm
    ccdc_input = l5to8m.select(bands_ccdc)

    # Format CCDC data values for optimization
    # Multiply by 10,000 and convert to int16 for ease of storage and computation
    # Copy image properties, including system properties
    ccdc_input = ccdc_input.map(format_ls_imagery)

    # Apply the CCDC algorithm using user-selected parameters
    ccdc = ee.Algorithms.TemporalSegmentation.Ccdc(
        collection =  ccdc_input,
        breakpointBands = breakpointBands,
        tmaskBands = tmaskBands,
        dateFormat = dateFormat,
        maxIterations = maxIterations # Changed to match global asset, maybe runs a lot faster?
#         lambda = lam -- lambda is reserved word in python, default is 20 which is what we are using anyway
    ) \
    .set({
      'start_year': start_year,
      'end_year': end_year,
      'start_doy': start_doy,
      'end_doy': end_doy,
      'breakpointBands': breakpointBands,
      'tmaskBands': tmaskBands,
      'dateFormat': dateFormat,
      'maxIterations': maxIterations,
      'maxCloudCoverLand': maxCloudCoverLand,
      'lambda': lam
    })

    # EXPORT -----

    # Export day of year composite
    task_doy = ee.batch.Export.image.toAsset(
        image = doys,
        description = 'ccdc_doys_' + tile_name + '_' + str(start_year) + '_' + str(end_year),
        assetId = asset_dir + 'CCDC_C2_SR_' + clipping_method + '_'+ str(start_year) + '_' + str(end_year) + '_' + version + \
        '_DOYs/CCDC_C2_' + str(start_year) + '_' + str(end_year) + '_' + tile_name + '_' + version + '_DOYs',
        region = clipper.geometry(),
        crs = 'EPSG:' + epsg,
        crsTransform = [30.0, 0.0, 15.0, 0.0, -30.0, 15],
        maxPixels = 1e13,
    )

    # Export CCDC model fits
    task_ccdc = ee.batch.Export.image.toAsset(
        image = ccdc,
        description = 'ccdc_' + tile_name + '_' + str(start_year) + '_' + str(end_year),
        assetId = asset_dir + 'CCDC_C2_SR_' + clipping_method + '_'+ str(start_year) + '_' + str(end_year) + '_' + version + \
        '/CCDC_C2_' + str(start_year) + '_' + str(end_year) + '_' + tile_name + '_' + version,
        region = clipper.geometry(),
        crs = 'EPSG:' + epsg,
        crsTransform = [30.0, 0.0, 15.0, 0.0, -30.0, 15],
        maxPixels = 1e13,
        pyramidingPolicy = {'.default': 'sample'}
    )

    task_doy.start()
    task_ccdc.start()

    print(tile_name, ': tasks submitted')

    return ccdc
