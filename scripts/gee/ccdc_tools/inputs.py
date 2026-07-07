"""
@license
Copyright 2019 Boston University

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Utility functions for getting inputs for CCDC.
"""

import ee
from . import dates as date_utils
from . import ccdc as ccdc_utils


# =========================
# 1. FUNCTIONS ==============
# =========================

def getLandsat(options):
    """
    Get Landsat images for a specific region.
    Possible bands and indices: BLUE, GREEN, RED, NIR, SWIR1, SWIR2, NDVI, NBR,
    EVI, EVI2, BRIGHTNESS, GREENNESS, WETNESS.

    @param {dict} options  Parameter dictionary containing:
        collection  {Number}      Landsat collection to use (1 or 2)
        start       {String}      First date to filter images
        end         {String}      Last date to filter images
        startDOY    {String}      First day of year to filter images
        endDOY      {String}      Last day of year to filter images
        region      {ee.Geometry} Region to filter the collection
        targetBands {list}        Bands and indices to return
        useMask     {Bool}        Whether to apply cloud/shadow masking
        sensors     {dict}        Dictionary with sensors to retrieve
    @returns {ee.ImageCollection} Landsat image collection filtered with the specified parameters
    """
    collection = (options and options.get('collection')) or 2
    start = (options and options.get('start')) or '1980-01-01'
    end = (options and options.get('end')) or '2023-01-01'
    start_doy = (options and options.get('startDOY')) or 1
    end_doy = (options and options.get('endDOY')) or 366
    region = (options and options.get('region')) or None
    target_bands = (options and options.get('targetBands')) or [
        'BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP',
        'NBR', 'NDFI', 'NDVI', 'GV', 'NPV', 'Shade', 'Soil',
        'EVI', 'EVI2', 'BRIGHTNESS', 'GREENNESS', 'WETNESS'
    ]
    use_mask = (options and options.get('useMask')) or True
    sensors = (options and options.get('sensors')) or {'l4': True, 'l5': True, 'l7': True, 'l8': True}

    if use_mask == 'No':
        use_mask = False

    # Define collection to use and select band names and functions accordingly
    if collection == 1:
        print('Landsat collection 1 has been deprecated')
    elif collection == 2:
        collection8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterDate(start, end)
        collection7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterDate(start, end)
        collection5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterDate(start, end)
        collection4 = ee.ImageCollection('LANDSAT/LT04/C02/T1_L2').filterDate(start, end)

        if use_mask:
            collection8 = collection8.map(prepareL8Col2)
            collection7 = collection7.map(prepareL4L5L7Col2)
            collection5 = collection5.map(prepareL4L5L7Col2)
            collection4 = collection4.map(prepareL4L5L7Col2)
        else:
            band_list_l8 = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'SR_B10']
            name_list_l8 = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP']
            band_list_l457 = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6']
            name_list_l457 = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP']
            collection8 = collection8.map(lambda i: i.select(band_list_l8).rename(name_list_l8))
            collection7 = collection7.map(lambda i: i.select(band_list_l457).rename(name_list_l457))
            collection4 = collection4.map(lambda i: i.select(band_list_l457).rename(name_list_l457))
            collection5 = collection5.map(lambda i: i.select(band_list_l457).rename(name_list_l457))

    # Merge all collections, compute indices and filter if requested
    col = collection4.merge(collection5).merge(collection7).merge(collection8)
    if region:
        col = col.filterBounds(region)

    indices = doIndices(col).select(target_bands)

    if not sensors.get('l5', True):
        indices = indices.filterMetadata('SATELLITE', 'not_equals', 'LANDSAT_5')
    if not sensors.get('l4', True):
        indices = indices.filterMetadata('SATELLITE', 'not_equals', 'LANDSAT_4')
    if not sensors.get('l7', True):
        indices = indices.filterMetadata('SATELLITE', 'not_equals', 'LANDSAT_7')
    if not sensors.get('l8', True):
        indices = indices.filterMetadata('SATELLITE', 'not_equals', 'LANDSAT_8')

    indices = indices.filter(ee.Filter.dayOfYear(start_doy, end_doy))
    return ee.ImageCollection(indices)


def doIndices(collection):
    """
    Calculate spectral indices for all bands in collection.
    @param {ee.ImageCollection} collection  Landsat image collection
    @returns {ee.ImageCollection} Landsat image with spectral indices
    """
    def add_indices(image):
        ndvi = calcNDVI(image)
        nbr = calcNBR(image)
        evi = calcEVI(image)
        evi2 = calcEVI2(image)
        tc = tcTrans(image)
        # NDFI function requires surface reflectance bands only
        bands = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']
        ndfi = calcNDFI(image.select(bands))
        return image.addBands([ndvi, nbr, evi, evi2, tc, ndfi])
    return collection.map(add_indices)


def getS2_old(roi):
    """
    Get Sentinel-2 surface reflectance data. Taken directly from GEE examples repo.
    @param {ee.Geometry} roi  Target study area to filter data
    @returns {ee.ImageCollection} Sentinel-2 SR and spectral indices
    """
    # Sentinel-2 Level 1C data. Bands B7, B8, B8A and B10 from this dataset are needed
    # as input to CDI and the cloud mask function.
    s2 = ee.ImageCollection('COPERNICUS/S2')
    # Cloud probability dataset. The probability band is used in the cloud mask function.
    s2c = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
    # Sentinel-2 surface reflectance data for the composite.
    s2_sr = ee.ImageCollection('COPERNICUS/S2_SR')

    # Dates over which to create a median composite.
    start = ee.Date('2016-03-01')
    end = ee.Date('2021-09-01')

    # S2 L1C for Cloud Displacement Index (CDI) bands.
    s2 = s2.filterDate(start, end).select(['B7', 'B8', 'B8A', 'B10'])
    # S2Cloudless for the cloud probability band.
    s2c = s2c.filterDate(start, end)
    # S2 L2A for surface reflectance bands.
    s2_sr = s2_sr.filterDate(start, end).select(['B2', 'B3', 'B4', 'B5', 'B8', 'B11', 'B12'])
    if roi:
        s2 = s2.filterBounds(roi)
        s2c = s2c.filterBounds(roi)
        s2_sr = s2_sr.filterBounds(roi)

    def index_join(collection_a, collection_b, property_name):
        joined = ee.ImageCollection(ee.Join.saveFirst(property_name).apply(**{
            'primary': collection_a,
            'secondary': collection_b,
            'condition': ee.Filter.equals(**{
                'leftField': 'system:index',
                'rightField': 'system:index'
            })
        }))
        # Merge the bands of the joined image.
        return joined.map(lambda image: image.addBands(ee.Image(image.get(property_name))))

    def mask_image(image):
        # Compute the cloud displacement index from the L1C bands.
        cdi = ee.Algorithms.Sentinel2.CDI(image)
        s2c_prob = image.select('probability')
        cirrus = image.select('B10').multiply(0.0001)
        # Assume low-to-mid atmospheric clouds to be pixels where probability
        # is greater than 65%, and CDI is less than -0.5.
        is_cloud = s2c_prob.gt(65).And(cdi.lt(-0.5)).Or(cirrus.gt(0.01))
        is_cloud = is_cloud.focal_min(3).focal_max(16)
        is_cloud = is_cloud.reproject(**{'crs': cdi.projection(), 'scale': 20})
        shadow_azimuth = ee.Number(90).subtract(
            ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE')))
        is_cloud = is_cloud.directionalDistanceTransform(shadow_azimuth, 50)
        is_cloud = is_cloud.reproject(**{'crs': cdi.projection(), 'scale': 100})
        is_cloud = is_cloud.select('distance').mask()
        return (image.select('B2', 'B3', 'B4', 'B8', 'B11', 'B12')
                     .rename(['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'])
                     .divide(10000).updateMask(is_cloud.Not())
                     .set('system:time_start',
                          ee.Image(image.get('l1c')).get('system:time_start')))

    # Join the cloud probability dataset to surface reflectance.
    with_cloud_probability = index_join(s2_sr, s2c, 'cloud_probability')
    # Join the L1C data to get the bands needed for CDI.
    with_s2l1c = index_join(with_cloud_probability, s2, 'l1c')
    # Map the cloud masking function over the joined collection.
    return doIndices(ee.ImageCollection(with_s2l1c.map(mask_image)))


def calcNDVI(image):
    """
    Calculate NDVI for an image.
    @param {ee.Image} image  Landsat image with NIR and RED bands
    @returns {ee.Image} NDVI image
    """
    ndvi = ee.Image(image).normalizedDifference(['NIR', 'RED']).rename('NDVI')
    return ndvi


def calcNBR(image):
    """
    Calculate NBR for an image.
    @param {ee.Image} image  Landsat image with NIR and SWIR2 bands
    @returns {ee.Image} NBR image
    """
    nbr = ee.Image(image).normalizedDifference(['NIR', 'SWIR2']).rename('NBR')
    return nbr


def calcNDFI(image):
    """
    Calculate NDFI using endmembers from Souza et al., 2005.
    @param {ee.Image} image  Surface reflectance image with 6 bands (i.e. not thermal)
    @returns {ee.Image} NDFI transform
    """
    gv = [.0500, .0900, .0400, .6100, .3000, .1000]
    shade = [0, 0, 0, 0, 0, 0]
    npv = [.1400, .1700, .2200, .3000, .5500, .3000]
    soil = [.2000, .3000, .3400, .5800, .6000, .5800]
    cloud = [.9000, .9600, .8000, .7800, .7200, .6500]
    cf = .1  # Not parameterized
    cf_threshold = ee.Image.constant(cf)
    unmix_image = (ee.Image(image).unmix([gv, shade, npv, soil, cloud], True, True)
                                  .rename(['band_0', 'band_1', 'band_2', 'band_3', 'band_4']))
    new_image = ee.Image(image).addBands(unmix_image)
    mask = new_image.select('band_4').lt(cf_threshold)
    ndfi = ee.Image(unmix_image).expression(
        '((GV / (1 - SHADE)) - (NPV + SOIL)) / ((GV / (1 - SHADE)) + NPV + SOIL)',
        {
            'GV': ee.Image(unmix_image).select('band_0'),
            'SHADE': ee.Image(unmix_image).select('band_1'),
            'NPV': ee.Image(unmix_image).select('band_2'),
            'SOIL': ee.Image(unmix_image).select('band_3')
        })
    return (ee.Image(new_image)
              .addBands(ee.Image(ndfi).rename(['NDFI']))
              .select(['band_0', 'band_1', 'band_2', 'band_3', 'NDFI'])
              .rename(['GV', 'Shade', 'NPV', 'Soil', 'NDFI'])
              .updateMask(mask))


def calcEVI(image):
    """
    Calculate EVI for an image.
    @param {ee.Image} image  Landsat image with NIR, RED, and BLUE bands
    @returns {ee.Image} EVI transform
    """
    evi = ee.Image(image).expression(
        'float(2.5*(((B4) - (B3)) / ((B4) + (6 * (B3)) - (7.5 * (B1)) + 1)))',
        {
            'B4': ee.Image(image).select(['NIR']),
            'B3': ee.Image(image).select(['RED']),
            'B1': ee.Image(image).select(['BLUE'])
        }).rename('EVI')
    return evi


def calcEVI2(image):
    """
    Calculate EVI2 for an image.
    @param {ee.Image} image  Landsat image with NIR and RED
    @returns {ee.Image} EVI2 transform
    """
    evi2 = ee.Image(image).expression(
        'float(2.5*(((B4) - (B3)) / ((B4) + (2.4 * (B3)) + 1)))',
        {
            'B4': image.select('NIR'),
            'B3': image.select('RED')
        })
    return evi2.rename('EVI2')


def tcTrans(image):
    """
    Tassel Cap coefficients from Crist 1985.
    @param {ee.Image} image  Landsat image with BLUE, GREEN, RED, NIR, SWIR1, and SWIR2
    @returns {ee.Image} 3-band image with Brightness, Greenness, and Wetness
    """
    brightness = image.expression(
        '(L1 * B1) + (L2 * B2) + (L3 * B3) + (L4 * B4) + (L5 * B5) + (L6 * B6)',
        {
            'L1': image.select('BLUE'), 'B1': 0.2043,
            'L2': image.select('GREEN'), 'B2': 0.4158,
            'L3': image.select('RED'), 'B3': 0.5524,
            'L4': image.select('NIR'), 'B4': 0.5741,
            'L5': image.select('SWIR1'), 'B5': 0.3124,
            'L6': image.select('SWIR2'), 'B6': 0.2303
        })
    greenness = image.expression(
        '(L1 * B1) + (L2 * B2) + (L3 * B3) + (L4 * B4) + (L5 * B5) + (L6 * B6)',
        {
            'L1': image.select('BLUE'), 'B1': -0.1603,
            'L2': image.select('GREEN'), 'B2': -0.2819,
            'L3': image.select('RED'), 'B3': -0.4934,
            'L4': image.select('NIR'), 'B4': 0.7940,
            'L5': image.select('SWIR1'), 'B5': -0.0002,
            'L6': image.select('SWIR2'), 'B6': -0.1446
        })
    wetness = image.expression(
        '(L1 * B1) + (L2 * B2) + (L3 * B3) + (L4 * B4) + (L5 * B5) + (L6 * B6)',
        {
            'L1': image.select('BLUE'), 'B1': 0.0315,
            'L2': image.select('GREEN'), 'B2': 0.2021,
            'L3': image.select('RED'), 'B3': 0.3102,
            'L4': image.select('NIR'), 'B4': 0.1594,
            'L5': image.select('SWIR1'), 'B5': -0.6806,
            'L6': image.select('SWIR2'), 'B6': -0.6109
        })
    bright = ee.Image(brightness).rename('BRIGHTNESS')
    green = ee.Image(greenness).rename('GREENNESS')
    wet = ee.Image(wetness).rename('WETNESS')
    tasseled_cap = ee.Image([bright, green, wet])
    return tasseled_cap


def makeLatGrid(min_y, max_y, min_x, max_x, size):
    """
    Create a grid with features corresponding to latitudinal strips.
    @param {Number} min_y  Minimum latitude coordinate
    @param {Number} max_y  Maximum latitude coordinate
    @param {Number} min_x  Minimum longitude coordinate
    @param {Number} max_x  Maximum longitude coordinate
    @param {Number} size   Size of features in units of latitudinal degrees
    @returns {ee.FeatureCollection} Grid of features along latitudinal lines
    """
    y_seq = ee.List.sequence(min_y, max_y, size)
    num_feats = y_seq.length().subtract(2)
    feats = ee.List.sequence(0, num_feats).map(lambda num: _lat_strip(num, y_seq, min_x, max_x))
    return ee.FeatureCollection(feats)


def _lat_strip(num, y_seq, min_x, max_x):
    num = ee.Number(num)
    num2 = num.add(1)
    y1 = ee.Number(y_seq.get(num))
    y2 = ee.Number(y_seq.get(num2))
    return ee.Feature(ee.Geometry.Polygon([[max_x, y2], [min_x, y2], [min_x, y1], [max_x, y1]]))


def makeLonGrid(min_y, max_y, min_x, max_x, size):
    """
    Create a grid with features corresponding to longitudinal strips.
    @param {Number} min_y  Minimum latitude coordinate
    @param {Number} max_y  Maximum latitude coordinate
    @param {Number} min_x  Minimum longitude coordinate
    @param {Number} max_x  Maximum longitude coordinate
    @param {Number} size   Size of features in units of longitudinal degrees
    @returns {ee.FeatureCollection} Grid of features along longitudinal lines
    """
    x_seq = ee.List.sequence(min_x, max_x, size)
    num_feats = x_seq.length().subtract(2)
    feats = ee.List.sequence(0, num_feats).map(lambda num: _lon_strip(num, x_seq, min_y, max_y))
    return ee.FeatureCollection(feats)


def _lon_strip(num, x_seq, min_y, max_y):
    num = ee.Number(num)
    num2 = num.add(1)
    x1 = ee.Number(x_seq.get(num))
    x2 = ee.Number(x_seq.get(num2))
    return ee.Feature(ee.Geometry.Polygon([[x2, max_y], [x1, max_y], [x1, min_y], [x2, min_y]]))


def makeLonLatGrid(min_y, max_y, min_x, max_x, size):
    """
    Create a grid with features along both latitudinal and longitudinal lines.
    @param {Number} min_y  Minimum latitude coordinate
    @param {Number} max_y  Maximum latitude coordinate
    @param {Number} min_x  Minimum longitude coordinate
    @param {Number} max_x  Maximum longitude coordinate
    @param {Number} size   Size of features in units of degrees
    @returns {ee.FeatureCollection} Grid of features
    """
    x_seq = ee.List.sequence(min_x, max_x, size)
    y_seq = ee.List.sequence(min_y, max_y, size)
    num_feats_y = y_seq.length().subtract(2)
    num_feats_x = x_seq.length().subtract(2)

    def make_row(y):
        y = ee.Number(y)
        y2 = y.add(1)
        y1_val = ee.Number(y_seq.get(y))
        y2_val = ee.Number(y_seq.get(y2))
        def make_col(x):
            x = ee.Number(x)
            x2 = x.add(1)
            x1_val = ee.Number(x_seq.get(x))
            x2_val = ee.Number(x_seq.get(x2))
            return ee.Feature(ee.Geometry.Polygon(
                [[x2_val, y2_val], [x1_val, y2_val], [x1_val, y1_val], [x2_val, y1_val]]))
        return ee.List.sequence(0, num_feats_x).map(make_col)

    feats = ee.List.sequence(0, num_feats_y).map(make_row)
    return ee.FeatureCollection(feats.flatten())


def makeAutoGrid(geo, size):
    """
    Create a grid with features overlaying the bounding box of a geometry.
    @param {ee.Geometry} geo   Geometry to use as spec for grid
    @param {Number}      size  Size of features in units of degrees
    @returns {ee.FeatureCollection} Grid of features
    """
    coord_list = ee.List(geo.coordinates().get(0))
    lon_list = coord_list.map(lambda c: ee.List(c).flatten().get(0))
    lat_list = coord_list.map(lambda c: ee.List(c).flatten().get(1))

    min_y = lat_list.reduce(ee.Reducer.min())
    max_y = lat_list.reduce(ee.Reducer.max())
    min_x = lon_list.reduce(ee.Reducer.min())
    max_x = lon_list.reduce(ee.Reducer.max())

    x_seq = ee.List.sequence(min_x, max_x, size)
    y_seq = ee.List.sequence(min_y, max_y, size)
    num_feats_y = y_seq.length().subtract(2)
    num_feats_x = x_seq.length().subtract(2)

    def make_row(y):
        y = ee.Number(y)
        y2 = y.add(1)
        y1_val = ee.Number(y_seq.get(y))
        y2_val = ee.Number(y_seq.get(y2))
        def make_col(x):
            x = ee.Number(x)
            x2 = x.add(1)
            x1_val = ee.Number(x_seq.get(x))
            x2_val = ee.Number(x_seq.get(x2))
            return ee.Feature(ee.Geometry.Polygon(
                [[x2_val, y2_val], [x1_val, y2_val], [x1_val, y1_val], [x2_val, y1_val]]))
        return ee.List.sequence(0, num_feats_x).map(make_col)

    feats = ee.List.sequence(0, num_feats_y).map(make_row)
    return ee.FeatureCollection(feats.flatten())


def getAncillary():
    """
    Get ancillary data for training and classification.
    @returns {ee.Image} Multi-band image containing ancillary layers
    """
    srtm = ee.Image('USGS/SRTMGL1_003').rename('ELEVATION')
    alos = ee.Image('JAXA/ALOS/AW3D30/V2_2').select(0).rename('ELEVATION')
    dem_image = ee.ImageCollection([alos, srtm]).mosaic()

    slope = ee.Terrain.slope(dem_image).rename('DEM_SLOPE')
    aspect = ee.Terrain.aspect(dem_image).rename('ASPECT')
    bio = ee.Image('WORLDCLIM/V1/BIO').select(['bio01', 'bio12']).rename(['TEMPERATURE', 'RAINFALL'])
    water = ee.Image('JRC/GSW1_1/GlobalSurfaceWater').select('occurrence').rename('WATER_OCCURRENCE')
    pop = (ee.ImageCollection('WorldPop/GP/100m/pop')
              .filterMetadata('year', 'equals', 2000)
              .mosaic()
              .rename('POPULATION'))
    hansen = (ee.Image('UMD/hansen/global_forest_change_2018_v1_6')
                .select('treecover2000')
                .rename('TREE_COVER'))
    night_lights = (ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG')
                       .filter(ee.Filter.date('2019-01-01', '2019-12-31'))
                       .select('avg_rad')
                       .mosaic()
                       .rename('NIGHT_LIGHTS'))

    return ee.Image.cat([dem_image, slope, aspect, bio, pop, water, hansen, night_lights]).unmask()


def maskS2clouds(image):
    """
    Mask Sentinel-2 imagery using QA band.
    @param {ee.Image} image  Sentinel-2 reflectance image
    @returns {ee.Image} Masked image
    """
    qa = image.select('QA60')
    # Bits 10 and 11 are clouds and cirrus, respectively.
    cloud_bit_mask = 1 << 10
    cirrus_bit_mask = 1 << 11
    # Both flags should be set to zero, indicating clear conditions.
    mask = (qa.bitwiseAnd(cloud_bit_mask).eq(0)
              .And(qa.bitwiseAnd(cirrus_bit_mask).eq(0)))
    return (image.updateMask(mask)
                 .select('B2', 'B3', 'B4', 'B8', 'B11', 'B12')
                 .rename(['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'])
                 .divide(10000)
                 .copyProperties(image, ['system:time_start']))


def getS2(roi):
    """
    Get Sentinel-2 surface reflectance data. Taken directly from GEE examples repo.
    @param {ee.Geometry} roi  Target study area to filter data
    @returns {ee.ImageCollection} Sentinel-2 SR and spectral indices
    """
    s2_sr = ee.ImageCollection('COPERNICUS/S2_SR')
    if roi:
        s2_sr = s2_sr.filterBounds(roi)
    s2_sr = s2_sr.map(maskS2clouds)
    return doIndices(s2_sr)


def getS1(focal_size=3, kernel_type='circle'):
    """
    Get Sentinel 1 data.
    @param {number} focal_size   Window size for focal mean (1 means no averaging)
    @param {string} kernel_type  Kernel type for focal mean
    @returns {ee.ImageCollection} Sentinel 1 collection with VH, VV, and ratio bands
    """
    def process_img(img):
        fmean = img.select('V.').add(30).focal_mean(focal_size, kernel_type)
        ratio = fmean.select('VH').divide(fmean.select('VV')).rename('ratio').multiply(30)
        angle = img.select('angle')
        return img.select().addBands(fmean).addBands(ratio).addBands(angle)

    data = (ee.ImageCollection('COPERNICUS/S1_GRD')
              .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
              .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
              .filter(ee.Filter.eq('instrumentMode', 'IW'))
              .select(['V.', 'angle'])
              .map(process_img))
    return data


def prepareL4L5(image):
    """
    Prepare Landsat 4 and 5 Collection 1 with strict filtering of noisy pixels.
    COLLECTION 1 HAS BEEN DEPRECATED.
    @param {ee.Image} image  Landsat SR image with pixel_qa band
    @returns {ee.Image} Landsat image with masked noisy pixels
    """
    band_list = ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6']
    name_list = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP']
    scaling = [10000, 10000, 10000, 10000, 10000, 10000, 1000]
    scaled = ee.Image(image).select(band_list).rename(name_list).divide(ee.Image.constant(scaling))
    valid_qa = [66, 130, 68, 132]
    mask1 = ee.Image(image).select(['pixel_qa']).remap(valid_qa, ee.List.repeat(1, len(valid_qa)), 0)
    mask2 = image.select('radsat_qa').eq(0)
    mask3 = image.select(band_list).reduce(ee.Reducer.min()).gt(0)
    mask4 = image.select('sr_atmos_opacity').unmask().lt(300)
    return ee.Image(image).addBands(scaled).updateMask(mask1.And(mask2).And(mask3).And(mask4))


def prepareL7(image):
    """
    Prepare Landsat 7 Collection 1 with strict filtering of noisy pixels.
    COLLECTION 1 HAS BEEN DEPRECATED.
    @param {ee.Image} image  Landsat SR image with pixel_qa band
    @returns {ee.Image} Landsat image with masked noisy pixels
    """
    band_list = ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6']
    name_list = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP']
    scaling = [10000, 10000, 10000, 10000, 10000, 10000, 1000]
    scaled = ee.Image(image).select(band_list).rename(name_list).divide(ee.Image.constant(scaling))
    valid_qa = [66, 130, 68, 132]
    mask1 = ee.Image(image).select(['pixel_qa']).remap(valid_qa, ee.List.repeat(1, len(valid_qa)), 0)
    mask2 = image.select('radsat_qa').eq(0)
    mask3 = image.select(band_list).reduce(ee.Reducer.min()).gt(0)
    mask4 = image.select('sr_atmos_opacity').unmask().lt(300)
    mask5 = ee.Image(image).mask().reduce(ee.Reducer.min()).focal_min(2.5)
    return ee.Image(image).addBands(scaled).updateMask(
        mask1.And(mask2).And(mask3).And(mask4).And(mask5))


def prepareL8(image):
    """
    Prepare Landsat 8 Collection 1 with strict filtering of noisy pixels.
    COLLECTION 1 HAS BEEN DEPRECATED.
    @param {ee.Image} image  Landsat SR image with pixel_qa band
    @returns {ee.Image} Landsat image with masked noisy pixels
    """
    band_list = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10']
    name_list = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP']
    scaling = [10000, 10000, 10000, 10000, 10000, 10000, 1000]
    valid_toa = [66, 68, 72, 80, 96, 100, 130, 132, 136, 144, 160, 164]
    valid_qa = [322, 386, 324, 388, 836, 900]
    scaled = ee.Image(image).select(band_list).rename(name_list).divide(ee.Image.constant(scaling))
    mask1 = ee.Image(image).select(['pixel_qa']).remap(valid_qa, ee.List.repeat(1, len(valid_qa)), 0)
    mask2 = image.select('radsat_qa').eq(0)
    mask3 = image.select(band_list).reduce(ee.Reducer.min()).gt(0)
    mask4 = ee.Image(image).select(['sr_aerosol']).remap(
        valid_toa, ee.List.repeat(1, len(valid_toa)), 0)
    return ee.Image(image).addBands(scaled).updateMask(mask1.And(mask2).And(mask3).And(mask4))


def prepareL4L5L7Col2(image):
    """
    Prepare Collection 2 Landsat 4, 5, and 7 with strict filtering of noisy pixels.
    @param {ee.Image} image  Landsat SR image with pixel_qa band
    @returns {ee.Image} Landsat image with masked noisy pixels
    """
    band_list = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6']
    name_list = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP']
    sub_band = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']
    optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
    thermal_band = image.select('ST_B6').multiply(0.00341802).add(149.0)
    scaled = optical_bands.addBands(thermal_band, None, True).select(band_list).rename(name_list)
    valid_qa = [5440, 5504]
    mask1 = ee.Image(image).select(['QA_PIXEL']).remap(valid_qa, ee.List.repeat(1, len(valid_qa)), 0)
    mask2 = image.select('QA_RADSAT').eq(0)
    mask3 = scaled.select(sub_band).reduce(ee.Reducer.min()).gt(0)
    mask4 = scaled.select(sub_band).reduce(ee.Reducer.max()).lt(1)
    # Mask hazy pixels using AOD threshold
    mask5 = image.select('SR_ATMOS_OPACITY').unmask(-1).lt(300)
    return ee.Image(image).addBands(scaled).updateMask(
        mask1.And(mask2).And(mask3).And(mask4).And(mask5))


def prepareL8Col2(image):
    """
    Prepare Collection 2 Landsat 8 with strict filtering of noisy pixels.
    @param {ee.Image} image  Landsat SR image with pixel_qa band
    @returns {ee.Image} Landsat image with masked noisy pixels
    """
    band_list = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10']
    name_list = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP']
    sub_band = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']
    optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
    thermal_band = image.select('ST_B10').multiply(0.00341802).add(149.0)
    scaled = optical_bands.addBands(thermal_band, None, True).select(band_list).rename(name_list)
    valid_toa = [2, 4, 32, 66, 68, 96, 100, 130, 132, 160, 164]
    valid_qa = [21824, 21888]
    mask1 = ee.Image(image).select(['QA_PIXEL']).remap(valid_qa, ee.List.repeat(1, len(valid_qa)), 0)
    mask2 = image.select('QA_RADSAT').eq(0)
    mask3 = scaled.select(sub_band).reduce(ee.Reducer.min()).gt(0)
    mask4 = scaled.select(sub_band).reduce(ee.Reducer.max()).lt(1)
    mask5 = ee.Image(image).select(['SR_QA_AEROSOL']).remap(
        valid_toa, ee.List.repeat(1, len(valid_toa)), 0)
    return ee.Image(image).addBands(scaled).updateMask(
        mask1.And(mask2).And(mask3).And(mask4).And(mask5))


def generateCollection(geom, start_date, end_date, collection=2):
    """
    Generate and combine filtered collections of Landsat 4, 5, 7 and 8.
    Simpler and faster than getLandsat.
    @param {ee.Geometry} geom        Geometry used to filter the collection
    @param {String}      start_date  Initial date to filter the collection
    @param {String}      end_date    Final date to filter the collection
    @param {Integer}     collection  Landsat collection to use (1 or 2)
    @returns {ee.ImageCollection} Filtered Landsat collection
    """
    if collection == 1:
        print('Collection 1 has been deprecated')
    elif collection == 2:
        filtered_l8 = (ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
                          .filter('WRS_ROW < 122')
                          .filterBounds(geom)
                          .map(prepareL8Col2))
        filtered_l7 = (ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
                          .filter('WRS_ROW < 122')
                          .filterBounds(geom)
                          .map(prepareL4L5L7Col2))
        filtered_l4 = (ee.ImageCollection('LANDSAT/LT04/C02/T1_L2')
                          .filter('WRS_ROW < 122')
                          .filterBounds(geom)
                          .map(prepareL4L5L7Col2))
        filtered_l5 = (ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
                          .filter('WRS_ROW < 122')
                          .filterBounds(geom)
                          .map(prepareL4L5L7Col2))

    merged_collections = (ee.ImageCollection(filtered_l8)
                            .merge(filtered_l7)
                            .merge(filtered_l5)
                            .merge(filtered_l4)
                            .filterDate(start_date, end_date))
    return merged_collections


def makeCcdImage(metadata_filter=None, segs=None, number_of_segments=None,
                  band_names=None, input_features=None, version=None):
    """
    Make a ccd image from the most recent known global run.
    @param {String} metadata_filter     Which ccdc run prefix to use
    @param {List}   segs                List with the segment names
    @param {Number} number_of_segments  Max number of segments to retrieve from the CCDC results
    @param {List}   band_names          List with the band names to use
    @param {List}   input_features      List with the CCDC features to extract
    @param {String} version             CCDC dataset version string
    @returns {ee.Image} Filtered CCDC results in 'long' format
    """
    metadata_filter = metadata_filter or 'z'
    number_of_segments = number_of_segments or 6
    band_names = band_names or ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP']
    segs = segs or ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']
    input_features = input_features or ['INTP', 'SLP', 'PHASE', 'AMPLITUDE', 'RMSE']
    version = version or 'v2'

    ccdc_collection = ee.ImageCollection('projects/CCDC/' + version)

    # Get CCDC coefficients
    ccdc_collection_filtered = ccdc_collection.filterMetadata(
        'system:index', 'starts_with', metadata_filter)

    # CCDC mosaic image
    ccdc = ccdc_collection_filtered.mosaic()

    # Turn array image into image
    return ee.Image(ccdc_utils.buildCcdImage(ccdc, number_of_segments, band_names))
