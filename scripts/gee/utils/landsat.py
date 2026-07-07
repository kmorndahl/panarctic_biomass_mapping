import ee
from . import landsat_proc_mod as lsat_proc

# FUNCTION: l4to8_c2_indices_xcal
# USE: Calculate a subset of spectral indices for which we have calculated cross-calibration parameters
# PARAMETERS:
#  img = pre-processed landsat image
# NOTES: map over image collection
# AUTHOR: Patrick Burns

def l4to8_c2_indices_xcal(img):

    # Normalized Difference Vegetation Index (NDVI)
    # Ref: Tucker, 1979
    NDVI = img.expression('((NIR - red) / (NIR + red))', {
    'NIR': img.select('NIR'),
    'red': img.select('red')
    }).rename('NDVI')

    # Enhanced Vegetation Index 2 band (EVI2b)
    # Ref: Jiang et al., 2008
    EVI2b = img.expression('(2.5 * (NIR - red) / (NIR + 2.5 * red + 1))', {
    'NIR': img.select('NIR'),
    'red': img.select('red'),
    'blue': img.select('blue')
    }).rename('EVI2b')

    # Normalized Burn Ratio (NBR)
    NBR = img.expression('((NIR - SWIR2) / (NIR + SWIR2))', {
    'NIR': img.select('NIR'),
    'SWIR2': img.select('SWIR2')
    }).rename('NBR')

    # Normalized Difference Moisture Index (NDMI)
    NDMI = img.expression('((NIR - SWIR1) / (NIR + SWIR1))', {
    'NIR': img.select('NIR'),
    'SWIR1': img.select('SWIR1')
    }).rename('NDMI')

    # Normalized Difference Water Index (NDWI)
    NDWI = img.expression('((green - NIR) / (green + NIR))', {
    'green': img.select('green'),
    'NIR': img.select('NIR')
    }).rename('NDWI')

    return ee.Image(img.addBands(NDVI) \
                       .addBands(EVI2b) \
                       .addBands(NBR) \
                       .addBands(NDMI) \
                       .addBands(NDWI))

# FUNCTION: getCal
# USE: Grabs appropriate cross-calibration parameters from dictionaries specified in 'applyCrossCalibration' functions
# PARAMETERS:
#  image = pre-processed landsat image
#  bandName = name of band for which to apply cross-calibration parameters, specified within 'applyCrossCalibration' scripts
#  factorNumber = specifies parameter number i.e. B0 (intercept), B1 (1st order polynomial), B2 (2nd order polynomial), B3 (3rd order polynomial)
# AUTHOR: Matt Macander

def getCal(factors, bandName, factorNumber):
    return ee.List(factors.get(bandName)).getNumber(factorNumber)

# FUNCTION: applyCal
# USE: Does the math to apply Landsat cross-calibrations
# PARAMETERS:
#  image = pre-processed landsat image
#  factors = cross-calibration parameters, pulled from dictionary within 'applyCrossCalibration' scripts
#  bandName = name of band for which to apply cross-calibration parameters, specified within 'applyCrossCalibration' scripts
# NOTES: map over image collection
# AUTHOR: Matt Macander

def applyCal(image, factors, bandName):

    # Get intercept parameter and create constant image
    # Get first order polynomial and multiply
    # Get second order polynomial and multiply
    # Get third order polynomial and multiply
    return ee.Image(getCal(factors, bandName, 0)).add(image.select(bandName).multiply(getCal(factors, bandName, 1))) \
                                                 .add(image.select(bandName).pow(2).multiply(getCal(factors, bandName, 2))) \
                                                 .add(image.select(bandName).pow(3).multiply(getCal(factors, bandName, 3))) \
                                                 .rename(bandName)

# FUNCTION: applyCrossCalibrationLS5
# USE: Apply cross-calibration to Landsat 5
# PARAMETERS:
#  image = pre-processed landsat image
# NOTES: map over image collection
# AUTHOR: Matt Macander
# TODO:
#  Pull cross-calibration parameters from table asset rather than hard-coding here
#  OR create Google Earth Engine objects with calibration parameters, create one applyCrossCalibration function, pull in appropriate dictionary based on user-input parameter

def applyCrossCalibrationLS5(image):

    # Dictionary of cross-calibration parameters
    factors = ee.Dictionary({
      'blue': [-0.0009, 0.9356, 0.1343, 0.0000],
      'green': [0.0013, 0.851, 0.3442, 0.0000],
      'red': [-0.0056, 1.0005, 0.0000, 0.0000],
      'NIR': [0.0019, 0.9854, 0.0000, 0.0000],
      'SWIR1': [-0.0041, 1.0687, -0.3977, 0.5442],
      'SWIR2': [0.0011, 0.9769, 0.0486, 0.0000],
      'EVI2b': [0.0047, 1.0232, 0.0000, 0.0000],
      'NBR': [0.0003, 1.0147, -0.0243, 0.0000],
      'NDMI': [0.0039, 1.0062, 0.0000, 0.0000],
      'NDVI': [0.0097, 1.0368, 0.0000, 0.0000],
      'NDWI': [-0.034, 1.033, 0.0519, 0.0000]
    })

    # Apply cross-calibration parameters
    cal = applyCal(image, factors, 'blue') \
         .addBands(applyCal(image, factors, 'green')) \
         .addBands(applyCal(image, factors, 'red')) \
         .addBands(applyCal(image, factors, 'NIR')) \
         .addBands(applyCal(image, factors, 'SWIR1')) \
         .addBands(applyCal(image, factors, 'SWIR2')) \
         .addBands(applyCal(image, factors, 'EVI2b')) \
         .addBands(applyCal(image, factors, 'NBR')) \
         .addBands(applyCal(image, factors, 'NDMI')) \
         .addBands(applyCal(image, factors, 'NDVI')) \
         .addBands(applyCal(image, factors, 'NDWI'))

    return image.addBands(cal, None, True)

# FUNCTION: applyCrossCalibrationLS8
# USE: Apply cross-calibration to Landsat 8/9
# PARAMETERS:
#  image = pre-processed landsat image
# NOTES: map over image collection
# AUTHOR: Matt Macander
# LAST UPDATE:
# TODO:
#  Pull cross-calibration parameters from table asset rather than hard-coding here

def applyCrossCalibrationLS8(image):

    # Dictionary of cross-calibration parameters
    factors = ee.Dictionary({
      'blue':  [0.0098, 1.0142, 0.0000, 0.0000],
      'green': [0.004, 1.0254, 0.0000, 0.0000],
      'red':   [0.0026, 1.1088, -0.5078, 0.826],
      'NIR':   [0.0233, 0.8115, 0.2393, 0.0000],
      'SWIR1': [0.0075, 0.9884, 0.0000, 0.0000],
      'SWIR2': [0.0059, 0.8746, 0.7391, -1.1305],
      'EVI2b': [-0.0025, 0.7878, 0.4396, -0.2729],
      'NBR': [-0.0098, 0.9703, 0.0631, 0.0000],
      'NDMI': [-0.0233, 0.978, 0.0818, 0.0000],
      'NDVI': [-0.0119, 0.792, 0.2272, 0.0000],
      'NDWI': [0.0269, 1.0057, 0.3434, 0.4678]
    })

    # Apply cross-calibration parameters
    cal = applyCal(image, factors, 'blue') \
          .addBands(applyCal(image, factors, 'green')) \
          .addBands(applyCal(image, factors, 'red')) \
          .addBands(applyCal(image, factors, 'NIR')) \
          .addBands(applyCal(image, factors, 'SWIR1')) \
          .addBands(applyCal(image, factors, 'SWIR2')) \
          .addBands(applyCal(image, factors, 'EVI2b')) \
          .addBands(applyCal(image, factors, 'NBR')) \
          .addBands(applyCal(image, factors, 'NDMI')) \
          .addBands(applyCal(image, factors, 'NDVI')) \
          .addBands(applyCal(image, factors, 'NDWI'))

    return image.addBands(cal, None, True)

# FUNCTION: permuteRadiometricError
# USE: Permute Landsat data with radiometric calibration error based on published values from the Landsat satellites
#  LS5 and LS7 cited by Logan: https://www.sciencedirect.com/science/article/pii/S0034425712000338#t0055
#  LS8 source by same author: https://www.mdpi.com/2072-4292/6/12/12275
#  LS8 additional source: https://www.mdpi.com/2072-4292/7/1/600
#  LS5: +-7%, LS7: +-5%, LS8: +-3%
# PARAMETERS:
#  image = pre-processed landsat image
# NOTES: map over image collection
# AUTHOR: Kathleen Orndahl
# LAST UPDATE:

def permuteRadiometricError(seed, bands):

    def wrap(image):

        # Get the satellite (spacecraft) ID
        sat_id = ee.String(image.get('SPACECRAFT_ID'))

        # Generate image using random normal distribution, select single random number from it
        # Because the same radiometric calibration coefficient is applied across all Landsat data,
        # we want an image of all one random number (rather than the number varying by pixel)
        # Do we need to .clamp(-1, 1) bc percentages are calibration standards the sensors need to abide by???
        # https://www.asprs.org/wp-content/uploads/pers/1997journal/jul/1997_jul_853-858.pdf
        rand = ee.Number(ee.Image.random(seed=seed, distribution='normal') \
                                 .sample(ee.Geometry.Point(0, 0)) \
                                 .first() \
                                 .get('random'))

        # Create permutations based on satellite ID
        perm = ee.Image.constant(ee.Algorithms.If(sat_id.equals('LANDSAT_8'), \
                                 ee.Number(rand.multiply(0.03)), \
                                 ee.Algorithms.If(sat_id.equals('LANDSAT_7'), \
                                                  ee.Number(rand.multiply(0.05)), \
                                                  ee.Number(rand.multiply(0.07)))))

        image_bands = image.select(bands)
        image_perm = image_bands.add(image_bands.multiply(perm))
        image_final = image.addBands(srcImg=image_perm, overwrite=True)

        return image_final

    return wrap

# FUNCTION: permuteXcalError
# USE: Permute Landsat data with cross-calibration error based on post calibration RMSE from Arctic specific cross-calibration
# PARAMETERS:
#  image = pre-processed landsat image
# NOTES: map over image collection
# AUTHOR: Kathleen Orndahl
# LAST UPDATE:

def permuteXcalError(seed):

    def wrap(image):

        # Get the satellite (spacecraft) ID
        sat_id = ee.String(image.get('SPACECRAFT_ID'))

        # Generate image using random normal distribution, select single random number from it
        # Because the same cross-calibration coefficient is applied across all Landsat data,
        # we want an image of all one random number (rather than the number varying by pixel)
        # Add to seed to make sure we are permuting with a different random number than the radiometric calibration permutation
        rand = ee.Number(ee.Image.random(seed=ee.Number(seed).add(42), distribution='normal') \
                                 .sample(ee.Geometry.Point(0, 0)) \
                                 .first() \
                                 .get('random'))

        # Assign xcal_rmse lookup column
        num = ee.Algorithms.If(sat_id.equals('LANDSAT_8'),  \
                               2, \
                               ee.Algorithms.If(sat_id.equals('LANDSAT_7'), \
                                                1, \
                                                0))

        # Landsat 5 / Landsat 7 / Landsat 8 cross-calibration error
        xcal_rmse = ee.Dictionary({ \
          'blue':  [0.007, 0, 0.009], \
          'green': [0.007, 0, 0.007], \
          'red':   [0.007, 0, 0.008], \
          'NIR':   [0.012, 0, 0.015], \
          'SWIR1': [0.012, 0, 0.013], \
          'SWIR2': [0.009, 0, 0.008], \
          'EVI2b': [0.017, 0, 0.023], \
          'NBR': [0.03, 0, 0.032], \
          'NDMI': [0.027, 0, 0.029], \
          'NDVI': [0.028, 0, 0.033], \
          'NDWI': [0.027, 0, 0.034] \
        })

        # Apply cross-calibration permutations
        image_perm = image \
        .addBands(srcImg=image.select('blue').add(rand.multiply(ee.List(xcal_rmse.get('blue')).getNumber(num))), overwrite=True) \
        .addBands(srcImg=image.select('green').add(rand.multiply(ee.List(xcal_rmse.get('green')).getNumber(num))), overwrite=True) \
        .addBands(srcImg=image.select('red').add(rand.multiply(ee.List(xcal_rmse.get('red')).getNumber(num))), overwrite=True) \
        .addBands(srcImg=image.select('NIR').add(rand.multiply(ee.List(xcal_rmse.get('NIR')).getNumber(num))), overwrite=True) \
        .addBands(srcImg=image.select('SWIR1').add(rand.multiply(ee.List(xcal_rmse.get('SWIR1')).getNumber(num))), overwrite=True) \
        .addBands(srcImg=image.select('SWIR2').add(rand.multiply(ee.List(xcal_rmse.get('SWIR2')).getNumber(num))), overwrite=True) \
        .addBands(srcImg=image.select('EVI2b').add(rand.multiply(ee.List(xcal_rmse.get('EVI2b')).getNumber(num))), overwrite=True) \
        .addBands(srcImg=image.select('NBR').add(rand.multiply(ee.List(xcal_rmse.get('NBR')).getNumber(num))), overwrite=True) \
        .addBands(srcImg=image.select('NDMI').add(rand.multiply(ee.List(xcal_rmse.get('NDMI')).getNumber(num))), overwrite=True) \
        .addBands(srcImg=image.select('NDVI').add(rand.multiply(ee.List(xcal_rmse.get('NDVI')).getNumber(num))), overwrite=True) \
        .addBands(srcImg=image.select('NDWI').add(rand.multiply(ee.List(xcal_rmse.get('NDWI')).getNumber(num))), overwrite=True)

        # Alternative to copyProperties: https://groups.google.com/g/google-earth-engine-developers/c/kjPV1YGcE9w/m/EftqZQXMAAAJ
        return image.select().addBands(image_perm)

    return wrap

# FUNCTION: add_doy
# USE: Adds day of year (doy) band to image
# PARAMETERS:
#  img = pre-processed landsat image
# NOTES: map over image collection
# AUTHOR: Matt Macander
# LAST UPDATE:

def add_doy(img):

    mask = img.select(0).mask()
    doy = ee.Image(img.date().getRelative('day','year').add(1)).uint16().rename('doy').updateMask(mask)
    return img.addBands(doy)

# FUNCTION: maskWater_jrc_ann
# USE: mask water based on the JRC Yearly Water Classification History, v1.4
# PARAMETERS:
#  mask_seasonal = whether or not to mask out seasonal, choose True or False
#  buffer_dist = distance over which to buffer the water mask, 0 = no buffer
# NOTES: map over image collection
# AUTHOR: Matt Macander
# LAST UPDATE:
# TODO:

def maskWater_jrc_ann(mask_seasonal, buffer_dist):

    def wrap(img):

        # Get the year the image was acquired
        year = ee.Number(img.date().get('year'))
        year = ee.Algorithms.If(ee.Number(year).gt(2021), 2021, year); # No 2022 JRC data, use 2021 for masking 2022 imagery

        # Get the JRC water classification for that year
        jrc_ann = ee.ImageCollection("JRC/GSW1_4/YearlyHistory").filter(ee.Filter.eq('year', year)).max();

        # Unmask
        jrc_unmask = jrc_ann.unmask()

        # Calculate mask based on user-input on whether to mask seasonal water or not
        # 0 = no data, 1 = no water, 2 = seasonal water, 3 = permanent water
        jrc_m = ee.Image(ee.Algorithms.If(mask_seasonal, jrc_unmask.lte(1), jrc_unmask.lte(2)))

        # Perform buffer, if specified
        jrc_m_b = ee.Image(ee.Algorithms.If(ee.Number(buffer_dist).eq(0),
                                            jrc_m,
                                            jrc_m.focalMin({'radius': buffer_dist, 'kernelType': 'square', 'units': 'meters'})))

        return img.updateMask(jrc_m_b)

    return wrap

# FUNCTION: maskWater_jrc_occurrence
# USE: mask water based on occurrence percentages (% of years with water present) from the JRC Global Surface Water Mapping Layers, v1.4
# PARAMETERS:
#  occurrence_percent = occurrence percentage (% of years with water present) threshold to be considered water, lower % = more water classifed, 0% is equivalent to the 'max_extent' band from the JRC dataset
#  buffer_dist = distance over which to buffer the water mask, 0 = no buffer
# NOTES: map over image collection
# AUTHOR: Kathleen Orndahl
# LAST UPDATE:
# TODO:
def maskWater_jrc_occurrence(occurrence_percent, buffer_dist):

    # https:#gis.stackexchange.com/questions/325259/mask-water-bodies-in-gee

    def wrap(img):

        # Get the JRC occurrence layer
        # Water = integer value indicating occurrence %, non-water or unmapped = masked
        # Apply threshold percent
        # Self mask to force all valid water occurrence percents to '1'
        # Unmask so that masked areas become zeros
        # Note that masked areas include oceans, so these become zeros
        occurrence = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select('occurrence').gte(occurrence_percent).selfMask().unmask(0)

        # Get Open Street Map Water Classification layer
        # https:#gee-community-catalog.Org/projects/osm_water/
        # Grab ocean only (1)
        # Mask so that ocean = 1, not ocean = 0
        # NOTE: JRC water metadata no longer works for identifying ocean
        #  - JRC/GSW1_0/Metadata outdated and inaccurate for Greenland
        #  - JRC/GSW1_4/Metadata has observations even in the ocean
        osm_water = ee.ImageCollection("projects/sat-io/open-datasets/OSM_waterLayer").mosaic()
        osm_ocean = osm_water.eq(1).selfMask().unmask(0); # 1 = ocean, 0 = not ocean

        # Get the maximum value between the 'occurrence' and 'osm_ocean' layers
        # Land: occurrence = 0, osm_ocean = 0, stays 0
        # Rivers/lakes: occurrence = 1, osm_ocean = 0, stays 1
        # Ocean: occurrence = 0, osm_ocean = 1, becomes 1
        water_mask = occurrence.max(osm_ocean)

        # Reverse mask to remove water
        water_mask = water_mask.eq(0)

        # Perform buffer, if specified
        water_mask_b = ee.Image(ee.Algorithms.If(ee.Number(buffer_dist).eq(0),
                                                 water_mask,
                                                 water_mask.focalMin(radius = buffer_dist, kernelType = 'square', units = 'meters')))

        return img.updateMask(water_mask_b)

    return wrap

# FUNCTION: maskOcean
# USE: mask ocean
# PARAMETERS:
# NOTES: map over image collection
# AUTHOR: Kathleen Orndahl
# LAST UPDATE:
# TODO:
def maskOcean(img):

    ocean_mask = ee.Image('projects/arctic-biomass-mapping/assets/ROIs/ocean_mask')

    return img.updateMask(ocean_mask)

# FUNCTION: maskOceanIce
# USE: mask ocean and Greenland icesheet
# PARAMETERS:
# NOTES: map over image collection
# AUTHOR: Kathleen Orndahl
# LAST UPDATE:
# TODO:
def maskOceanIce(img):

    ocean_ice_mask = ee.Image('projects/arctic-biomass-mapping/assets/ROIs/ocean_ice_mask_v20240206')

    return img.updateMask(ocean_ice_mask)
