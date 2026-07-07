import ee 

# Landsat Collection 2 Surface Reflectance Image Collections
l4 = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2") # July 1982 to December 1993
l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2") # March 1984 to May 2012
l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2") # April 1999 to present
l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2") # February 2013 to present
l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2") # October 2021 to present

# FUNCTION: l4to8_c2_rename_bands
# USE: rename Landsat 4,5,7,8 bands
# AUTHOR: Patrick Burns

def l4to8_c2_rename_bands(img):
    
    # Get the satellite (spacecraft) ID
    sat_id = ee.String(img.get('SPACECRAFT_ID'))

    # L4, L5, and L7 have the same bands. L8 is slightly different. Use if logic to rename
    new_band_names = ee.Algorithms.If(sat_id.compareTo('LANDSAT_8'),
    # new band names for L4, L5, and L7
    ["blue", "green", "red", "NIR", "SWIR1", "SWIR2",
    "SR_ATMOS_OPACITY", "SR_CLOUD_QA", "ST",
    "ST_ATRAN", "ST_CDIST", "ST_DRAD", "ST_EMIS", "ST_EMSD",
    "ST_QA", "ST_TRAD", "ST_URAD", "QA_PIXEL", "QA_RADSAT"],
    # new bands names for L8
    ["u_blue", "blue", "green", "red", "NIR", "SWIR1", "SWIR2",
    "SR_QA_AEROSOL", "ST", "ST_ATRAN", "ST_CDIST", "ST_DRAD", "ST_EMIS",
    "ST_EMSD", "ST_QA", "ST_TRAD", "ST_URAD", "QA_PIXEL", "QA_RADSAT"]
    )

    return img.rename(new_band_names)

# FUNCTION: l457_c2_scaleoff
# USE: Apply scale and offset values to convert Landsat 4,5,7 to reflectance values (0 to 1)
# AUTHOR: Patrick Burns
# NOTE: need to change band names first using function l4to8_c2_rename_bands

def l457_c2_scaleoff(img):
    
    # Save system properties
    sys_footprint = img.get("system:footprint")
    sys_time_start = img.get("system:time_start")
    sys_time_end = img.get("system:time_end")

    # These bands don't need to be rescaled
    noadj_bands = img.select(["SR_CLOUD_QA", "QA_PIXEL", "QA_RADSAT"])

    # Apply the reflectance scale and offset
    sr_bands_adj = img.select(["blue", "green", "red", "NIR", "SWIR1", "SWIR2"]).multiply(0.0000275).add(-0.2)

    # Apply the temp scale and offset
    st_band_adj = img.select(['ST']).multiply(0.00341802).add(149)

    # These bands only need to be scaled by 0.01
    s01_bands_adj = img.select(["ST_CDIST", "ST_QA"]).multiply(0.01)

    # These bands only need to be scaled by 0.001
    s001_bands_adj = img.select(["SR_ATMOS_OPACITY", "ST_DRAD", "ST_TRAD", "ST_URAD"]).multiply(0.001)

    # These bands only need to be scaled by 0.0001
    s0001_bands_adj = img.select(["ST_ATRAN", "ST_EMIS", "ST_EMSD"]).multiply(0.0001)

    return sr_bands_adj.addBands(st_band_adj) \
                       .addBands(s01_bands_adj) \
                       .addBands(s001_bands_adj) \
                       .addBands(s0001_bands_adj) \
                       .addBands(noadj_bands) \
                       .copyProperties(img) \
                       .set({'system:footprint': sys_footprint,
                             'system:time_start': sys_time_start,
                             'system:time_end': sys_time_end})

# FUNCTION: l8_c2_scaleoff
# USE: Apply scale and offset values to convert Landsat 8 to reflectance values (0 to 1)
# AUTHOR: Patrick Burns
# NOTE: need to change band names first using function l4to8_c2_rename

def l8_c2_scaleoff(img):
    
    # Save system properties
    sys_footprint = img.get("system:footprint")
    sys_time_start = img.get("system:time_start")
    sys_time_end = img.get("system:time_end")

    # These bands don't need to be rescaled
    noadj_bands = img.select(["SR_QA_AEROSOL", "QA_PIXEL", "QA_RADSAT"])

    # Apply the reflectance scale and offset
    sr_bands_adj = img.select(["u_blue", "blue", "green", "red", "NIR", "SWIR1", "SWIR2"]).multiply(0.0000275).add(-0.2)

    # Apply the temp scale and offset
    st_band_adj = img.select(['ST']).multiply(0.00341802).add(149)

    # These bands only need to be scaled by 0.01
    s01_bands_adj = img.select(["ST_CDIST", "ST_QA"]).multiply(0.01)

    # These bands only need to be scaled by 0.001
    s001_bands_adj = img.select(["ST_DRAD", "ST_TRAD", "ST_URAD"]).multiply(0.001)

    # These bands only need to be scaled by 0.0001
    s0001_bands_adj = img.select(["ST_ATRAN", "ST_EMIS", "ST_EMSD"]).multiply(0.0001)

    return sr_bands_adj.addBands(st_band_adj) \
                       .addBands(s01_bands_adj) \
                       .addBands(s001_bands_adj) \
                       .addBands(s0001_bands_adj) \
                       .addBands(noadj_bands) \
                       .copyProperties(img) \
                       .set({'system:footprint': sys_footprint,
                             'system:time_start': sys_time_start,
                             'system:time_end': sys_time_end})

# FUNCTION: l4to8_c2_scaleoff
# USE: Apply scale and offset values to convert Landsat 4,5,7,8 to reflectance values (0 to 1)
# AUTHOR: Patrick Burns

def l4to8_c2_scaleoff(img):
    
    # Get the satellite (spacecraft) ID
    sat_id = ee.String(img.get('SPACECRAFT_ID'))

    img_adj = ee.Algorithms.If(sat_id.compareTo('LANDSAT_8'),
    l457_c2_scaleoff(img),
    l8_c2_scaleoff(img))

    return(img_adj)

# FUNCTION: l4to8_c2_maskSR
# USE: Keep surface reflectance values greater than 0 and less than or equal to 1
# AUTHOR: Patrick Burns
# NOTE: blue and u_blue bands are excluded from the masking procedure since they tend to be noiser due to atmospheric effects

def l4to8_c2_maskSR(img):

    sat_id = ee.String(img.get('SPACECRAFT_ID'))
    sr_adj = ee.Image(ee.Algorithms.If(sat_id.compareTo('LANDSAT_8'),
                    # bands for for L4, L5, and L7 (don't worry about blue for the mask)
                       img.select(["blue", "SR_ATMOS_OPACITY", "SR_CLOUD_QA", "ST",
                                   "ST_ATRAN", "ST_CDIST", "ST_DRAD", "ST_EMIS", "ST_EMSD",
                                   "ST_QA", "ST_TRAD", "ST_URAD", "QA_PIXEL", "QA_RADSAT"]) \
                          .addBands(img.select(["green", "red", "NIR", "SWIR1", "SWIR2"]) \
                                       .updateMask(img.select(["green", "red", "NIR", "SWIR1", "SWIR2"]) \
                                                   .reduce(ee.Reducer.min()).gt(0) \
                                                   .And(img.select(["green", "red", "NIR", "SWIR1", "SWIR2"]) \
                                                           .reduce(ee.Reducer.max()).lte(1)))),
                    # bands for for L8 (don't worry about u_blue or blue for the mask)
                    img.select(["u_blue", "blue", "SR_QA_AEROSOL", "ST", "ST_ATRAN", "ST_CDIST", "ST_DRAD", "ST_EMIS",
                                "ST_EMSD", "ST_QA", "ST_TRAD", "ST_URAD", "QA_PIXEL", "QA_RADSAT"]) \
                       .addBands(img.select(["green", "red", "NIR", "SWIR1", "SWIR2"]) \
                                    .updateMask(img.select(["green", "red", "NIR", "SWIR1", "SWIR2"]) \
                                                   .reduce(ee.Reducer.min()).gt(0) \
                                                   .And(img.select(["green", "red", "NIR", "SWIR1", "SWIR2"]) \
                                                           .reduce(ee.Reducer.max()).lte(1))))
                    ))
    return sr_adj

# FUNCTION: l4to8_c2_indices
# USE: Calculate a set of spectral indices
# AUTHOR: Patrick Burns

def l4to8_c2_indices(img):
    
    # Normalized Difference Vegetation Index (NDVI)
    # Ref: Tucker, 1979
    NDVI = img.expression('((NIR - red) / (NIR + red))', {
    'NIR': img.select('NIR'),
    'red': img.select('red')
    }).rename('NDVI')

    # Kernel NDVI
    # Ref: Camps-Valls et al., 2021 (https:#doi.Org/10.1126/sciadv.abc7447)
    def kNDVI_calc(img):
        
        red = img.select('red')
        NIR = img.select('NIR')
        D2 = NIR.subtract(red).pow(2).select([0],['d2']); # Compute D2 and rename it to d2

        # --------------------------------------------------------------------------------------
        # Fix sigma parameter to a reasonable value for all pixels in the region or for your specific problem
        # We recommend values ranging between [0.1-0.3] when working with (normalized) reflectance data.
        sigma = 0.15

        # kDNVI
        return D2.divide(sigma).divide(sigma).divide(4.0).tanh().select([0], ['kNDVI'])

    kNDVI = kNDVI_calc(img)

    # Enhanced Vegetation Index (EVI)
    # Ref: Huete et al., 1997
    EVI = img.expression('(2.5 * (NIR - red) / (NIR + 6 * red - 7.5 * blue + 1))', {
    'NIR': img.select('NIR'),
    'red': img.select('red'),
    'blue': img.select('blue')
    }).rename('EVI')

    # Enhanced Vegetation Index 2 band (EVI2b)
    # Ref: Jiang et al., 2008
    EVI2b = img.expression('(2.5 * (NIR - red) / (NIR + 2.5 * red + 1))', {
    'NIR': img.select('NIR'),
    'red': img.select('red'),
    'blue': img.select('blue')
    }).rename('EVI2b')

    # Global Environmental Monitoring Index (GEMI)
    # Ref: Pinty and Verstraete, 1992
    GEMI = img.expression('((2 * ((NIR * NIR) - (red * red)) + 1.5 * NIR + 0.5 * red) / (NIR + red + 0.5)) * (1 - 0.25 * ((2 * ((NIR * NIR)     - (red * red)) + 1.5 * NIR + 0.5 * red) / (NIR + red + 0.5))) - (red - 0.125) / (1 - red)', {
    'red': img.select('red'),
    'NIR': img.select('NIR')
    }).rename('GEMI')

    # Tasseled cap transformation
    # Ref: Crist, 1985; Kennedy et al., 2010 (https:#doi.Org/10.1016/j.rse.2010.07.008)
    # Tasseled Cap Greenness (TCg)
    TCg = img.expression('-0.1603 * blue - 0.2819 * green - 0.4934 * red + 0.7940 * NIR - 0.0002 * SWIR1 - 0.1446 * SWIR2', {
    'blue': img.select('blue'),
    'green': img.select('green'),
    'red': img.select('red'),
    'NIR': img.select('NIR'),
    'SWIR1': img.select('SWIR1'),
    'SWIR2': img.select('SWIR2')
    }).rename('TCg')

    # Tasseled Cap Wetness (TCw)
    # Ref: Crist, 1985
    TCw = img.expression('0.0315 * blue + 0.2021 * green + 0.3102 * red + 0.1594 * NIR - 0.6806 * SWIR1 - 0.6109 * SWIR2', {
    'blue': img.select('blue'),
    'green': img.select('green'),
    'red': img.select('red'),
    'NIR': img.select('NIR'),
    'SWIR1': img.select('SWIR1'),
    'SWIR2': img.select('SWIR2')
    }).rename('TCw')

    # Tasseled Cap Brightness (TCb)
    TCb = img.expression('0.2043 * blue + 0.4158 * green + 0.5524 * red + 0.5741 * NIR + 0.3124 * SWIR1 + 0.2303 * SWIR2', {
    'blue': img.select('blue'),
    'green': img.select('green'),
    'red': img.select('red'),
    'NIR': img.select('NIR'),
    'SWIR1': img.select('SWIR1'),
    'SWIR2': img.select('SWIR2')
    }).rename('TCb')

    # Normalized difference fraction index (NDFI)
    # Ref: Souza et al., 2005
    def calc_NDFI(img):
        
        gv = [0.05, 0.09, 0.04, 0.61, 0.3, 0.1] # green veg
        shade = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        npv = [0.14, 0.17, 0.22, 0.3, 0.55, 0.3] # non-photosynthetic veg
        soil = [0.2, 0.3, 0.34, 0.58, 0.6, 0.58]
        cloud = [0.9, 0.96, 0.8, 0.78, 0.72, 0.65]
        unmixed = img.select(['blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2']).unmix({
          'endmembers': [gv, shade, npv, soil, cloud],
          'sumToOne': True,
          'nonNegative': True
        }) \
        .rename(['gv', 'shade', 'npv', 'soil', 'cloud'])

        return unmixed.expression(
          '((i.gv / (1 - i.shade)) - (i.npv + i.soil)) / ((i.gv / (1 - i.shade)) + i.npv + i.soil)',
          {'i': unmixed}
        )

    NDFI = calc_NDFI(img).rename('NDFI')

    # Spectral Variability Vegetation Index (SVVI)
    # Ref: Coulter et al., 2016
    def calc_SVVI(img):
        
        sd_1to7 = img.select(['blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2']).reduce(ee.Reducer.stdDev())
        sd_457 = img.select(['NIR', 'SWIR1', 'SWIR2']).reduce(ee.Reducer.stdDev())
        svvi = sd_1to7.subtract(sd_457)

        return svvi

    SVVI = calc_SVVI(img).rename('SVVI')

    # Soil Adjusted Vegetation Index (SAVI)
    # Ref: Huete, 1988
    SAVI = img.expression('(((NIR - red) / (NIR + red *0.5)) * (1 + 0.5))', {
    'NIR': img.select('NIR'),
    'red': img.select('red')
    }).rename('SAVI')

    # Modified Soil Adjusted Vegetation Index (MSAVI)
    MSAVI = img.expression('1/2 * ((2*(NIR+1)) - (((2*NIR)+1)**2 - 8*(NIR-red))**(1/2))', {
    'NIR': img.select('NIR'),
    'red': img.select('red'),
    }).rename('MSAVI')

    # Normalized Burn Ratio (NBR)
    NBR = img.expression('((NIR - SWIR2) / (NIR + SWIR2))', {
    'NIR': img.select('NIR'),
    'SWIR2': img.select('SWIR2')
    }).rename('NBR')

    # Normalized Burn Ratio 2 (NBR2)
    NBR2 = img.expression('((SWIR1 - SWIR2) / (SWIR1 + SWIR2))', {
    'SWIR1': img.select('SWIR1'),
    'SWIR2': img.select('SWIR2')
    }).rename('NBR2')

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

    # Normalized Difference Snow Index (NDSI)
    NDSI = img.expression('((green - SWIR1) / (green + SWIR1))', {
    'green': img.select('green'),
    'SWIR1': img.select('SWIR1')
    }).rename('NDSI')

    return ee.Image(img.addBands(NDVI)
                       .addBands(kNDVI)
                       .addBands(EVI)
                       .addBands(EVI2b)
                       .addBands(GEMI)
                       .addBands(TCg)
                       .addBands(TCw)
                       .addBands(TCb)
                       .addBands(NDFI)
                       .addBands(SVVI)
                       .addBands(SAVI)
                       .addBands(MSAVI)
                       .addBands(NBR)
                       .addBands(NBR2)
                       .addBands(NDMI)
                       .addBands(NDWI)
                       .addBands(NDSI))

# FUNCTION: l4to8_c2_indices_noblue
# USE: Calculate a set of spectral indices
# AUTHOR: Patrick Burns

def l4to8_c2_indices_noblue(img):
    
    # Normalized Difference Vegetation Index (NDVI)
    # Ref: Tucker, 1979
    NDVI = img.expression('((NIR - red) / (NIR + red))', {
    'NIR': img.select('NIR'),
    'red': img.select('red')
    }).rename('NDVI')

    # Kernel NDVI
    # Ref: Camps-Valls et al., 2021 (https:#doi.Org/10.1126/sciadv.abc7447)
    def kNDVI_calc(img):
        
        red = img.select('red')
        NIR = img.select('NIR')
        D2 = NIR.subtract(red).pow(2).select([0],['d2']); # Compute D2 and rename it to d2

        # --------------------------------------------------------------------------------------
        # Fix sigma parameter to a reasonable value for all pixels in the region or for your specific problem
        # We recommend values ranging between [0.1-0.3] when working with (normalized) reflectance data.
        sigma = 0.15

        # kDNVI
        return D2.divide(sigma).divide(sigma).divide(4.0).tanh().select([0], ['kNDVI'])

    kNDVI = kNDVI_calc(img)

    # Enhanced Vegetation Index (EVI)
    # Ref: Huete et al., 1997
    EVI = img.expression('(2.5 * (NIR - red) / (NIR + 6 * red - 7.5 * blue + 1))', {
    'NIR': img.select('NIR'),
    'red': img.select('red'),
    'blue': img.select('blue')
    }).rename('EVI')

    # Enhanced Vegetation Index 2 band (EVI2b)
    # Ref: Jiang et al., 2008
    EVI2b = img.expression('(2.5 * (NIR - red) / (NIR + 2.5 * red + 1))', {
    'NIR': img.select('NIR'),
    'red': img.select('red'),
    }).rename('EVI2b')

    # Global Environmental Monitoring Index (GEMI)
    # Ref: Pinty and Verstraete, 1992
    GEMI = img.expression('((2 * ((NIR * NIR) - (red * red)) + 1.5 * NIR + 0.5 * red) / (NIR + red + 0.5)) * (1 - 0.25 * ((2 * ((NIR * NIR)     - (red * red)) + 1.5 * NIR + 0.5 * red) / (NIR + red + 0.5))) - (red - 0.125) / (1 - red)', {
    'red': img.select('red'),
    'NIR': img.select('NIR')
    }).rename('GEMI')



    # Normalized difference fraction index (NDFI)
    # Ref: Souza et al., 2005
    def calc_NDFI(img):
        
        gv = [0.09, 0.04, 0.61, 0.3, 0.1] # green veg
        shade = [0.0, 0.0, 0.0, 0.0, 0.0]
        npv = [0.17, 0.22, 0.3, 0.55, 0.3] # non-photosynthetic veg
        soil = [0.3, 0.34, 0.58, 0.6, 0.58]
        cloud = [0.96, 0.8, 0.78, 0.72, 0.65]
        unmixed = img.select(['green', 'red', 'NIR', 'SWIR1', 'SWIR2']).unmix({
          'endmembers': [gv, shade, npv, soil, cloud],
          'sumToOne': True,
          'nonNegative': True
        }) \
        .rename(['gv', 'shade', 'npv', 'soil', 'cloud'])

        return unmixed.expression(
          '((i.gv / (1 - i.shade)) - (i.npv + i.soil)) / ((i.gv / (1 - i.shade)) + i.npv + i.soil)',
          {'i': unmixed}
        )

    NDFI = calc_NDFI(img).rename('NDFI')

    # Spectral Variability Vegetation Index (SVVI)
    # Ref: Coulter et al., 2016
    def calc_SVVI(img):

        sd_2to7 = img.select(['green', 'red', 'NIR', 'SWIR1', 'SWIR2']).reduce(ee.Reducer.stdDev())
        sd_457 = img.select(['NIR', 'SWIR1', 'SWIR2']).reduce(ee.Reducer.stdDev())
        svvi = sd_2to7.subtract(sd_457)

        return svvi

    SVVI = calc_SVVI(img).rename('SVVI')

    # Soil Adjusted Vegetation Index (SAVI)
    # Ref: Huete, 1988
    SAVI = img.expression('(((NIR - red) / (NIR + red *0.5)) * (1 + 0.5))', {
    'NIR': img.select('NIR'),
    'red': img.select('red')
    }).rename('SAVI')

    # Modified Soil Adjusted Vegetation Index (MSAVI)
    MSAVI = img.expression('1/2 * ((2*(NIR+1)) - (((2*NIR)+1)**2 - 8*(NIR-red))**(1/2))', {
    'NIR': img.select('NIR'),
    'red': img.select('red'),
    }).rename('MSAVI')

    # Normalized Burn Ratio (NBR)
    NBR = img.expression('((NIR - SWIR2) / (NIR + SWIR2))', {
    'NIR': img.select('NIR'),
    'SWIR2': img.select('SWIR2')
    }).rename('NBR')

    # Normalized Burn Ratio 2 (NBR2)
    NBR2 = img.expression('((SWIR1 - SWIR2) / (SWIR1 + SWIR2))', {
    'SWIR1': img.select('SWIR1'),
    'SWIR2': img.select('SWIR2')
    }).rename('NBR2')

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

    # Normalized Difference Snow Index (NDSI)
    NDSI = img.expression('((green - SWIR1) / (green + SWIR1))', {
    'green': img.select('green'),
    'SWIR1': img.select('SWIR1')
    }).rename('NDSI')

    return ee.Image(img.addBands(NDVI)
                       .addBands(kNDVI)
                       .addBands(EVI2b)
                       .addBands(GEMI)
                       .addBands(NDFI)
                       .addBands(SVVI)
                       .addBands(SAVI)
                       .addBands(MSAVI)
                       .addBands(NBR)
                       .addBands(NBR2)
                       .addBands(NDMI)
                       .addBands(NDWI)
                       .addBands(NDSI))

# FUNCTION: l4to8_c2_qa_maskClouds
# USE: mask clouds based on the QA_PIXEL band (derived from CFMask algorithm).
#      Also remove radiometrically-saturated pixels using QA_RADSAT band
# AUTHOR: Patrick Burns

def l4to8_c2_qa_maskClouds(img):

    # Calculate bitmask values
    fillBitMask = ee.Number(2).pow(0).int()
    dilatedCloudsBitMask = ee.Number(2).pow(1).int()
    cirrusBitMask = ee.Number(2).pow(2).int(); #only applicable for L8
    cloudsBitMask = ee.Number(2).pow(3).int()
    cloudShadowBitMask = ee.Number(2).pow(4).int()

    # For Landsat 8 USGS states: "pixels classified as high aerosol content are not recommended for use"
    # Use if logic to add the aersol mask for Landsat 8 only
    qp = img.select(['QA_PIXEL'])
    sat_id = ee.String(img.get('SPACECRAFT_ID'))
    qp_mask = ee.Image(ee.Algorithms.If(sat_id.compareTo('LANDSAT_8'),
             # mask for for L4, L5, and L7
             qp.bitwiseAnd(fillBitMask).eq(0) \
               .And(qp.bitwiseAnd(dilatedCloudsBitMask).eq(0)) \
               .And(qp.bitwiseAnd(cirrusBitMask).eq(0)) \
               .And(qp.bitwiseAnd(cloudsBitMask).eq(0)) \
               .And(qp.bitwiseAnd(cloudShadowBitMask).eq(0)),
             # mask for L8
             qp.bitwiseAnd(fillBitMask).eq(0) \
               .And(qp.bitwiseAnd(dilatedCloudsBitMask).eq(0)) \
               .And(qp.bitwiseAnd(cirrusBitMask).eq(0)) \
               .And(qp.bitwiseAnd(cloudsBitMask).eq(0)) \
               .And(qp.bitwiseAnd(cloudShadowBitMask).eq(0)) \
               .And(img.select('SR_QA_AEROSOL').gt(1)) \
               .And(img.select('SR_QA_AEROSOL').lte(164))
             ))

    # Also make a mask which excludes radiometrically saturated pixels
    radsat_mask = ee.Image(img.select('QA_RADSAT').eq(0))

    return img.updateMask(qp_mask.And(radsat_mask))

# FUNCTION: l4to8_c2_qa_maskWater
# USE: mask water based on the QA_PIXEL band (derived from CFMask algorithm).
# AUTHOR: Patrick Burns

def l4to8_c2_qa_maskWater(img):

    water_mask = ee.Image(img.select(['QA_PIXEL']).bitwiseAnd(ee.Number(2).pow(7).int()).eq(0))

    return img.updateMask(water_mask)

# FUNCTION: maskWater_jrc_ann
# USE: mask water based on the JRC Yearly Water Classification History, v1.3
# AUTHOR: Patrick Burns
# NOTE: v1.3 only goes through 2020

def maskWater_jrc_ann(img):

    # Get the year the image was acquired
    year = ee.Number(img.date().get('year'))

    # Get the JRC water classification for that year
    jrc_ann = ee.ImageCollection("JRC/GSW1_3/YearlyHistory").filter(ee.Filter.eq('year', year)).max()

    # Calculate the mask based on seasonal and permanent water and buffer the mask neighborhood using a radius of 30 m
    jrc_mask = jrc_ann.unmask().lte(1)
                        #.focalMin({radius: 30, kernelType: 'square', units: 'meters'})

    return img.updateMask(jrc_mask)

# FUNCTION: l4to8_c2_qa_maskSnow
# USE: mask snow based on the QA_PIXEL band (derived from CFMask algorithm).
# AUTHOR: Patrick Burns

def l4to8_c2_qa_maskSnow(img):

    water_mask = ee.Image(img.select(['QA_PIXEL']).bitwiseAnd(ee.Number(2).pow(5).int()).eq(0))

    return img.updateMask(water_mask)

# FUNCTION: l4to8_c2_topocorr_ic
# USE: compute the illumination condition of a scene using a DEM
# AUTHOR: Patrick Burns

def l4to8_c2_topocorr_ic(img, dem):

    # Smooth the DEM using a low-pass kernel
    boxcar = ee.Kernel.square({'radius': 3,
                                 'units': 'pixels',
                                 'normalize': True})

    dem_s = ee.Image(dem).convolve(boxcar)

    # Extract image metadata about solar position and convert to constant images
    SZ_rad = ee.Image.constant(ee.Number(90).subtract(ee.Number(img.get('SUN_ELEVATION')))) \
                     .multiply(3.14159265359) \
                     .divide(180) \
                     .clip(img.geometry().buffer(10000))
    SA_rad = ee.Image.constant(ee.Number(img.get('SUN_AZIMUTH')).multiply(3.14159265359).divide(180)).clip(img.geometry().buffer(10000))

    # Create terrain layers
    slp = ee.Terrain.slope(dem_s).clip(img.geometry().buffer(10000))
    slp_rad = ee.Terrain.slope(dem_s).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000))
    asp_rad = ee.Terrain.aspect(dem_s).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000))

    # Calculate the Illumination Condition (IC)
    # slope part of the illumination condition
    cosZ = SZ_rad.cos()
    cosS = slp_rad.cos()
    slope_illumination = cosS.expression("cosZ * cosS",
                                          {'cosZ': cosZ,
                                           'cosS': cosS.select('slope')})
    # aspect part of the illumination condition
    sinZ = SZ_rad.sin()
    sinS = slp_rad.sin()
    cosAziDiff = (SA_rad.subtract(asp_rad)).cos()
    aspect_illumination = sinZ.expression("sinZ * sinS * cosAziDiff",
                                           {'sinZ': sinZ,
                                            'sinS': sinS,
                                            'cosAziDiff': cosAziDiff})
    # full illumination condition (IC)
    ic = slope_illumination.add(aspect_illumination)

    # Add IC to original image
    img_ic = ee.Image(img.addBands(ic.rename('IC')) \
                         .addBands(cosZ.rename('cosZ')) \
                         .addBands(cosS.rename('cosS')) \
                         .addBands(slp.rename('slope')))
    
    return img_ic

# FUNCTION: l4to8_c2_topocorr_scsc
# USE: apply the sun canopy sensor topographic correction
# AUTHOR: Patrick Burns

def l4to8_c2_topocorr_scsc(imgcol, dem):

    def _tc(img):

        # Save system properties
        sys_footprint = img.get("system:footprint")
        sys_time_start = img.get("system:time_start")
        sys_time_end = img.get("system:time_end")

        # Add the illumination condition
        img_plus_ic = l4to8_c2_topocorr_ic(img, dem)

        # Build masks
        mask1 = img_plus_ic.select('NIR').gt(0)
        mask2 = img_plus_ic.select('slope').gte(5).And(img_plus_ic.select('IC').gt(0))

        img_plus_ic_mask2 = ee.Image(img_plus_ic.updateMask(mask2))

        # Specify Bands to topographically correct
        bandList = ['blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2']

        def apply_SCSccorr(band):
            
            method = 'SCSc'
            out = img_plus_ic_mask2.select('IC', band).reduceRegion({
              'reducer': ee.Reducer.linearFit(), # Compute coefficients: a(slope), b(offset), c(b/a)
              'geometry': ee.Geometry(img.geometry().buffer(-5000)), # trim off the outer edges of the image for linear relationship
              'scale': 30,
              'maxPixels': 1000000000
              })

            out_a = ee.Number(out.get('scale'))
            out_b = ee.Number(out.get('offset'))
            out_c = ee.Number(out.get('offset')).divide(ee.Number(out.get('scale')))

            #apply the SCSc correction
            SCSc_output = img_plus_ic_mask2.expression("((image * (cosB * cosZ + cvalue)) / (ic + cvalue))", {
              'image': img_plus_ic_mask2.select(band),
              'ic': img_plus_ic_mask2.select('IC'),
              'cosB': img_plus_ic_mask2.select('cosS'),
              'cosZ': img_plus_ic_mask2.select('cosZ'),
              'cvalue': out_c
              })

            return ee.Image(SCSc_output)

        # Apply the topo correction over each spectral band
        img_SCSccorr = ee.Image(bandList.map(apply_SCSccorr)).addBands(img_plus_ic.select('IC'))

        # Erode the illumination condition mask
        ic_mask = img_plus_ic.select('IC').gt(0).focalMin({radius: 30, kernelType: 'square', units: 'meters'})

        return img_SCSccorr.unmask(img_plus_ic.select(ee.List([bandList, "IC"]).flatten())) \
                           .updateMask(ic_mask) \
                           .addBands(img_plus_ic.select('u_.*|SR_.*|ST_.*|QA_.*')) \
                           .copyProperties(img) \
                           .set({'system:footprint': sys_footprint,
                                 'system:time_start': sys_time_start,
                                 'system:time_end': sys_time_end})

    return imgcol.map(_tc)

# FUNCTION: l4to8_c2_topocorr_se
# USE: apply the Statistical Empirical topographic correction
# AUTHOR: Patrick Burns

def l4to8_c2_topocorr_se(imgcol, dem):

    def _tc(img):

        # Save system properties
        sys_footprint = img.get("system:footprint")
        sys_time_start = img.get("system:time_start")
        sys_time_end = img.get("system:time_end")

        # Add the illumination condition
        img_plus_ic = l4to8_c2_topocorr_ic(img, dem)

        # Build masks
        mask1 = img_plus_ic.select('NIR').gt(0)
        mask2 = img_plus_ic.select('slope').gte(5).And(img_plus_ic.select('IC').gt(0))

        img_plus_ic_mask2 = ee.Image(img_plus_ic.updateMask(mask2))

        # Specify Bands to topographically correct
        bandList = ['blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2']

        def apply_SEcorr(band):
            
            method = 'SE'
            out = img_plus_ic_mask2.select('IC', band).reduceRegion({
            'reducer': ee.Reducer.linearFit(), # Compute coefficients: a(slope), b(offset), c(b/a)
            'geometry': ee.Geometry(img.geometry().buffer(-5000)),
            'scale': 30,
            'maxPixels': 1000000000
            }) 

            out_a = ee.Number(out.get('scale'))
            out_b = ee.Number(out.get('offset'))
            out_c = ee.Number(out.get('offset')).divide(ee.Number(out.get('scale')))

            out2 = img_plus_ic_mask2.select(band).reduceRegion({
            'reducer': ee.Reducer.mean(), # Compute mean value for each band
            'geometry': ee.Geometry(img.geometry().buffer(-5000)),
            'scale': 30,
            'maxPixels': 1000000000
            })

            out2mean = ee.Number(out2.get(band))

            #apply the SE correction
            se_output = img_plus_ic_mask2.expression("((image - (avalue * ic + bvalue) + bandmean))", {
              'image': img_plus_ic_mask2.select(band),
              'ic': img_plus_ic_mask2.select('IC'),
              'cosZ': img_plus_ic_mask2.select('cosZ'),
              'avalue': out_a,
              'bvalue': out_b,
              'bandmean': out2mean
            });#.rename(bandList+'corr_'+method)

            keyTopoCorrMeth = 'TopoCorr_Method'
            keya = band+'a'
            keyb = band+'b'
            keyc = band+'c'

            return ee.Image(se_output).setMulti({'TopoCorr_Meth': method, 'TopoCorr_a': out_a, 'TopoCorr_b': out_b, 'TopoCorr_c': out_c}) \
                                      .copyProperties(img_plus_ic_mask2)

        img_SEcorr = ee.Image(bandList.map(apply_SEcorr)).addBands(img_plus_ic.select('IC'))

        # Erode the illumination condition mask
        ic_mask = img_plus_ic.select('IC').gt(0).focalMin({radius: 30, kernelType: 'square', units: 'meters'})

        return img_SEcorr.unmask(img_plus_ic.select(ee.List([bandList, "IC"]).flatten())) \
                         .updateMask(ic_mask) \
                         .addBands(img_plus_ic.select('u_.*|SR_.*|ST_.*|QA_.*')) \
                         .copyProperties(img) \
                         .set({'system:footprint': sys_footprint,
                               'system:time_start': sys_time_start,
                               'system:time_end': sys_time_end})

    return imgcol.map(_tc)
