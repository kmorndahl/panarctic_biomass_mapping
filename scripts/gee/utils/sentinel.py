"""

DESCRIPTION: Functions for processing Sentinel-2 data

Merged from:
  - sentinel2_proc_mod.py (Patrick Burns, pb463@nau.edu, NAU)
  - 00_fun_sentinel.py (Melissa Rose, Kathleen Orndahl, NAU)

FUNCTION LIST:
  scale_refl_toFlt
  sel_rename_bands
  scale_refl_toUint16
  L1C_filter_col
  L2A_filter_col
  L2A_mask_cloudless_wshadow
  L1C_topo_corr
  add_spec_indices
  brdf_corr
  quality_S2_CloudlessComposite

"""

import ee


# =========================
# 1. REFLECTANCE SCALING =
# =========================

# FUNCTION: scale_refl_toFlt
# USE: rescale reflectance bands to floating point values, so reflectance values range from [0,1]
# NOTE: only keeps reflectance bands; works for L1C and L2A

def scale_refl_toFlt(img):
    date = img.date()
    return (ee.Image(img.select('B.*')
                       .divide(10000)
                       .copyProperties(img))
               .set('system:time_start', date.millis()))


# FUNCTION: sel_rename_bands
# USE: rename select S2 bands
# NOTE: works for L1C and L2A

def sel_rename_bands(img):
    orig_sel_names = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']
    new_sel_names  = ['Blue', 'Green', 'Red', 'RE1', 'RE2', 'RE3', 'NIR', 'RE4', 'SWIR1', 'SWIR2']
    return ee.Image(img.select(orig_sel_names)
                       .rename(new_sel_names))


# FUNCTION: scale_refl_toUint16
# USE: rescale reflectance bands to Uint16, so reflectance values range from [0,65535]
# NOTE: only keeps reflectance bands; works for L1C and L2A

def scale_refl_toUint16(img):
    return ee.Image(img.select('B.*')
                       .uint16()
                       .copyProperties(img))


# ==========================
# 2. COLLECTION FILTERING ==
# ==========================

# FUNCTION: L1C_filter_col
# USE: filter S2 image collection using a geometry (aoi), date range, and maximum cloud percentage (cloudy_perc_filter)

def L1C_filter_col(aoi, start_date, end_date, cloud_perc_filter):
    return (ee.ImageCollection('COPERNICUS/S2')
              .filterBounds(aoi)
              .filterDate(start_date, end_date)
              .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', cloud_perc_filter)))


# FUNCTION: L2A_filter_col
# USE: filter S2 image collection using a geometry (aoi), date range, and maximum cloud percentage (cloudy_perc_filter)

def L2A_filter_col(aoi, start_date, end_date, cloud_perc_filter):
    return (ee.ImageCollection('COPERNICUS/S2_SR')
              .filterBounds(aoi)
              .filterDate(start_date, end_date)
              .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', cloud_perc_filter)))


# ==============================
# 3. CLOUD / SHADOW MASKING ====
# ==============================

# FUNCTION: L2A_mask_cloudless_wshadow
# USE: apply S2 cloudless mask and remove shadows
# ref: https://developers.google.com/earth-engine/tutorials/community/sentinel-2-s2cloudless

def L2A_mask_cloudless_wshadow(S2_L2A_col, cloud_prob_thresh, NIR_dark_thresh, cloud_proj_dist, buffer):

    S2_cloudless_col = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')

    # Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
    S2_wCld = ee.ImageCollection(
        ee.Join.saveFirst('s2cloudless').apply({
            'primary':   S2_L2A_col,
            'secondary': S2_cloudless_col,
            'condition': ee.Filter.equals({
                'leftField':  'system:index',
                'rightField': 'system:index'
            })
        })
    )

    def add_cloud_bands(img):
        cld_prb = ee.Image(img.get('s2cloudless')).select('probability')
        is_cloud = cld_prb.gt(cloud_prob_thresh).rename('clouds')
        return img.addBands(ee.Image([cld_prb, is_cloud]))

    def add_shadow_bands(img):
        not_water = img.select('SCL').neq(6)
        SR_band_scale = 1e4
        dark_pixels = (img.select('B8')
                          .lt(NIR_dark_thresh * SR_band_scale)
                          .multiply(not_water)
                          .rename('dark_pixels'))
        shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')))
        cld_proj = (img.select('clouds')
                       .directionalDistanceTransform(shadow_azimuth, cloud_proj_dist * 10)
                       .reproject({'crs': img.select(0).projection(), 'scale': 100})
                       .select('distance')
                       .mask()
                       .rename('cloud_transform'))
        shadows = cld_proj.multiply(dark_pixels).rename('shadows')
        return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))

    def add_cld_shdw_mask(img):
        img_cloud = add_cloud_bands(img)
        img_cloud_shadow = add_shadow_bands(img_cloud)
        is_cld_shdw = (img_cloud_shadow.select('clouds')
                                       .add(img_cloud_shadow.select('shadows'))
                                       .gt(0))
        is_cld_shdw = (is_cld_shdw.focal_min(2)
                                   .focal_max(buffer * 2 / 20)
                                   .reproject({'crs': img.select([0]).projection(), 'scale': 20})
                                   .rename('cloudmask'))
        return img.addBands(is_cld_shdw)

    def apply_cld_shdw_mask(img):
        not_cld_shdw = img.select('cloudmask').Not()
        return img.select('B.*').updateMask(not_cld_shdw)

    return (S2_wCld.map(add_cld_shdw_mask)
                   .map(apply_cld_shdw_mask))


# ==============================
# 4. TOPOGRAPHIC CORRECTION ====
# ==============================

# FUNCTION: L1C_topo_corr
# USE: apply topographic correction to an L2A image
# NOTE: need to rename bands first using sel_rename_bands()

def L1C_topo_corr(img):
    """Apply Sun-Canopy-Sensor + C (SCSc) topographic correction to an image.

    Functions by Patrick Burns (pb463@nau.edu) and Matt Macander (mmacander@abrinc.com).
    """

    def illuminationCondition(img):
        SZ_rad = (ee.Image.constant(ee.Number(img.get('MEAN_SOLAR_ZENITH_ANGLE')))
                    .multiply(3.14159265359)
                    .divide(180)
                    .clip(img.geometry().buffer(10000)))
        SA_rad = (ee.Image.constant(
                      ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
                        .multiply(3.14159265359)
                        .divide(180))
                    .clip(img.geometry().buffer(10000)))
        dem = ee.Image('NASA/NASADEM_HGT/001')
        slp     = ee.Terrain.slope(dem).clip(img.geometry().buffer(10000))
        slp_rad = ee.Terrain.slope(dem).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000))
        asp_rad = ee.Terrain.aspect(dem).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000))
        cosZ = SZ_rad.cos()
        cosS = slp_rad.cos()
        slope_illumination = cosS.expression(
            'cosZ * cosS',
            {'cosZ': cosZ, 'cosS': cosS.select('slope')}
        )
        sinZ = SZ_rad.sin()
        sinS = slp_rad.sin()
        cosAziDiff = (SA_rad.subtract(asp_rad)).cos()
        aspect_illumination = sinZ.expression(
            'sinZ * sinS * cosAziDiff',
            {'sinZ': sinZ, 'sinS': sinS, 'cosAziDiff': cosAziDiff}
        )
        ic = slope_illumination.add(aspect_illumination)
        img_plus_ic = ee.Image(
            img.addBands(ic.rename('IC'))
               .addBands(cosZ.rename('cosZ'))
               .addBands(cosS.rename('cosS'))
               .addBands(slp.rename('slope'))
        )
        return img_plus_ic

    def illuminationCorrection(img):
        props = img.toDictionary()
        st    = img.get('system:time_start')
        img_plus_ic = img
        mask1 = img_plus_ic.select('NIR').gt(-0.1)
        mask2 = (img_plus_ic.select('slope').gte(5)
                            .And(img_plus_ic.select('IC').gte(0))
                            .And(img_plus_ic.select('NIR').gt(-0.1)))
        img_plus_ic_mask2 = ee.Image(img_plus_ic.updateMask(mask2))
        bandList = ['Blue', 'Green', 'Red', 'RE1', 'RE2', 'RE3', 'NIR', 'RE4', 'SWIR1', 'SWIR2']
        compositeBands  = img.bandNames()
        nonCorrectBands = img.select(compositeBands.removeAll(bandList))

        def apply_SCSccorr(band):
            out = (ee.Image(1)
                     .addBands(img_plus_ic_mask2.select('IC', band))
                     .reduceRegion(
                         reducer=ee.Reducer.linearRegression(2, 1),
                         geometry=ee.Geometry(img.geometry()),
                         scale=300,
                         bestEffort=True,
                         maxPixels=1e10
                     ))
            fit = out.combine({'coefficients': ee.Array([[1], [1]])}, False)
            out_a = ee.Array(fit.get('coefficients')).get([0, 0])
            out_b = ee.Array(fit.get('coefficients')).get([1, 0])
            out_c = out_a.divide(out_b)
            SCSc_output = img_plus_ic_mask2.expression(
                '((image * (cosB * cosZ + cvalue)) / (ic + cvalue))',
                {
                    'image':  img_plus_ic_mask2.select(band),
                    'ic':     img_plus_ic_mask2.select('IC'),
                    'cosB':   img_plus_ic_mask2.select('cosS'),
                    'cosZ':   img_plus_ic_mask2.select('cosZ'),
                    'cvalue': out_c
                }
            )
            return SCSc_output

        img_SCSccorr  = ee.Image([apply_SCSccorr(b) for b in bandList]).addBands(img_plus_ic.select('IC'))
        bandList_IC   = ee.List([bandList, 'IC']).flatten()
        img_SCSccorr  = img_SCSccorr.unmask(img_plus_ic.select(bandList_IC)).select(bandList)
        return ee.Image(img_SCSccorr).setMulti(props).set('system:time_start', st)

    img_plus_ic = illuminationCondition(img)
    return illuminationCorrection(img_plus_ic)


# =========================
# 5. SPECTRAL INDICES =====
# =========================

# FUNCTION: add_spec_indices
# USE: calculate spectral indices
# NOTE: use scaled reflectance values [0-1] from scale_refl_toFlt(); rename bands using sel_rename_bands(); works for L1C and L2A

def add_spec_indices(img):

    NDVI = ee.Image(img.normalizedDifference(['NIR', 'Red'])
                       .select([0], ['NDVI']))

    def kNDVI_calc(img):
        Red = img.select('Red')
        NIR = img.select('NIR')
        D2  = NIR.subtract(Red).pow(2).select([0], ['d2'])
        sigma = 0.15
        return D2.divide(sigma).divide(sigma).divide(4.0).tanh().select([0], ['kNDVI'])

    kNDVI = kNDVI_calc(img)

    def EVI_calc(img):
        Red  = img.select('Red')
        NIR  = img.select('NIR')
        Blue = img.select('Blue')
        num = NIR.subtract(Red).multiply(ee.Image(2.5))
        den = NIR.add(Red.multiply(6.0)).subtract(Blue.multiply(7.5)).add(1.0)
        return (num.divide(den)).select([0], ['EVI'])

    EVI = EVI_calc(img)

    RGVI  = img.normalizedDifference(['Red', 'Green']).select([0], ['RGVI'])
    NDMI  = img.normalizedDifference(['NIR', 'SWIR1']).select([0], ['NDMI'])
    NBR   = img.normalizedDifference(['NIR', 'SWIR2']).select([0], ['NBR'])
    CIRE  = (img.select('RE3')
               .divide(img.select('RE1'))
               .subtract(ee.Image.constant(1))
               .rename(['CIRE']))

    def SVVI_calc(img):
        term1 = ee.Image(img.select(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])).reduce(ee.Reducer.stdDev())
        term2 = ee.Image(img.select(['NIR', 'SWIR1', 'SWIR2'])).reduce(ee.Reducer.stdDev())
        return (term1.subtract(term2)).select([0], ['SVVI'])

    SVVI = SVVI_calc(img)

    return (img.addBands(NDVI)
               .addBands(kNDVI)
               .addBands(EVI)
               .addBands(RGVI)
               .addBands(NDMI)
               .addBands(NBR)
               .addBands(CIRE)
               .addBands(SVVI))


# =========================
# 6. BRDF CORRECTION ======
# =========================

# FUNCTION: brdf_corr
# USE: correct for BRDF
# ref: https://doi.org/10.3390/rs11070831

def brdf_corr(img):

    PI = ee.Number(3.14159265359)
    MAX_SATELLITE_ZENITH = 7.5
    MAX_DISTANCE = 1000000
    UPPER_LEFT  = 0
    LOWER_LEFT  = 1
    LOWER_RIGHT = 2
    UPPER_RIGHT = 3

    def value(lst, index):
        return ee.Number(lst.get(index))

    def where(condition, trueValue, falseValue):
        trueMasked  = trueValue.mask(condition)
        falseMasked = falseValue.mask(invertMask(condition))
        return trueMasked.unmask(falseMasked)

    def invertMask(mask):
        return mask.multiply(-1).add(1)

    def line_from_coords(coordinates, fromIndex, toIndex):
        return ee.Geometry.LineString(ee.List([
            coordinates.get(fromIndex),
            coordinates.get(toIndex)
        ]))

    def getsunAngles(date, footprint):
        jdp = date.getFraction('year')
        seconds_in_hour = 3600
        hourGMT = ee.Number(date.getRelative('second', 'day')).divide(seconds_in_hour)
        latRad  = ee.Image.pixelLonLat().select('latitude').multiply(PI.divide(180))
        longDeg = ee.Image.pixelLonLat().select('longitude')
        jdpr = jdp.multiply(PI).multiply(2)
        a = ee.List([0.000075, 0.001868, 0.032077, 0.014615, 0.040849])
        meanSolarTime = longDeg.divide(15.0).add(ee.Number(hourGMT))
        localSolarDiff1 = (value(a, 0)
                           .add(value(a, 1).multiply(jdpr.cos()))
                           .subtract(value(a, 2).multiply(jdpr.sin()))
                           .subtract(value(a, 3).multiply(jdpr.multiply(2).cos()))
                           .subtract(value(a, 4).multiply(jdpr.multiply(2).sin())))
        localSolarDiff2 = localSolarDiff1.multiply(12 * 60)
        localSolarDiff  = localSolarDiff2.divide(PI)
        trueSolarTime   = (meanSolarTime
                           .add(localSolarDiff.divide(60))
                           .subtract(12.0))
        ah = trueSolarTime.multiply(ee.Number(MAX_SATELLITE_ZENITH * 2).multiply(PI.divide(180)))
        b  = ee.List([0.006918, 0.399912, 0.070257, 0.006758, 0.000907, 0.002697, 0.001480])
        delta = (value(b, 0)
                 .subtract(value(b, 1).multiply(jdpr.cos()))
                 .add(value(b, 2).multiply(jdpr.sin()))
                 .subtract(value(b, 3).multiply(jdpr.multiply(2).cos()))
                 .add(value(b, 4).multiply(jdpr.multiply(2).sin()))
                 .subtract(value(b, 5).multiply(jdpr.multiply(3).cos()))
                 .add(value(b, 6).multiply(jdpr.multiply(3).sin())))
        cosSunZen = (latRad.sin().multiply(delta.sin())
                     .add(latRad.cos().multiply(ah.cos()).multiply(delta.cos())))
        sunZen = cosSunZen.acos()
        sinSunAzSW = ah.sin().multiply(delta.cos()).divide(sunZen.sin())
        sinSunAzSW = sinSunAzSW.clamp(-1.0, 1.0)
        cosSunAzSW = ((latRad.cos().multiply(-1).multiply(delta.sin())
                       .add(latRad.sin().multiply(delta.cos()).multiply(ah.cos())))
                      .divide(sunZen.sin()))
        sunAzSW = sinSunAzSW.asin()
        sunAzSW = where(cosSunAzSW.lte(0), sunAzSW.multiply(-1).add(PI), sunAzSW)
        sunAzSW = where(cosSunAzSW.gt(0).And(sinSunAzSW.lte(0)), sunAzSW.add(PI.multiply(2)), sunAzSW)
        sunAz = sunAzSW.add(PI)
        sunAz = where(sunAz.gt(PI.multiply(2)), sunAz.subtract(PI.multiply(2)), sunAz)
        footprint_polygon = ee.Geometry.Polygon(footprint)
        sunAz  = sunAz.clip(footprint_polygon).rename(['sunAz'])
        sunZen = sunZen.clip(footprint_polygon).rename(['sunZen'])
        return [sunAz, sunZen]

    def azimuth(footprint):
        def x(point): return ee.Number(ee.List(point).get(0))
        def y(point): return ee.Number(ee.List(point).get(1))
        upperCenter = line_from_coords(footprint, UPPER_LEFT, UPPER_RIGHT).centroid().coordinates()
        lowerCenter = line_from_coords(footprint, LOWER_LEFT, LOWER_RIGHT).centroid().coordinates()
        slope       = ((y(lowerCenter)).subtract(y(upperCenter))).divide((x(lowerCenter)).subtract(x(upperCenter)))
        slopePerp   = ee.Number(-1).divide(slope)
        azimuthLeft = ee.Image(PI.divide(2).subtract((slopePerp).atan()))
        return azimuthLeft.rename(['viewAz'])

    def zenith(footprint):
        leftLine      = line_from_coords(footprint, UPPER_LEFT, LOWER_LEFT)
        rightLine     = line_from_coords(footprint, UPPER_RIGHT, LOWER_RIGHT)
        leftDistance  = ee.FeatureCollection(leftLine).distance(MAX_DISTANCE)
        rightDistance = ee.FeatureCollection(rightLine).distance(MAX_DISTANCE)
        viewZenith    = (rightDistance.multiply(ee.Number(MAX_SATELLITE_ZENITH * 2))
                                      .divide(rightDistance.add(leftDistance))
                                      .subtract(ee.Number(MAX_SATELLITE_ZENITH))
                                      .clip(ee.Geometry.Polygon(footprint))
                                      .rename(['viewZen']))
        return viewZenith.multiply(PI.divide(180))

    def _apply(image, kvol, kvol0):
        blue  = _correct_band(image, 'Blue',  kvol, kvol0, f_iso=0.0774, f_geo=0.0079, f_vol=0.0372)
        green = _correct_band(image, 'Green', kvol, kvol0, f_iso=0.1306, f_geo=0.0178, f_vol=0.0580)
        red   = _correct_band(image, 'Red',   kvol, kvol0, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574)
        re1   = _correct_band(image, 'RE1',   kvol, kvol0, f_iso=0.2085, f_geo=0.0256, f_vol=0.0845)
        re2   = _correct_band(image, 'RE2',   kvol, kvol0, f_iso=0.2316, f_geo=0.0273, f_vol=0.1003)
        re3   = _correct_band(image, 'RE3',   kvol, kvol0, f_iso=0.2599, f_geo=0.0294, f_vol=0.1197)
        nir   = _correct_band(image, 'NIR',   kvol, kvol0, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535)
        re4   = _correct_band(image, 'RE4',   kvol, kvol0, f_iso=0.2907, f_geo=0.0410, f_vol=0.1611)
        swir1 = _correct_band(image, 'SWIR1', kvol, kvol0, f_iso=0.3430, f_geo=0.0453, f_vol=0.1154)
        swir2 = _correct_band(image, 'SWIR2', kvol, kvol0, f_iso=0.2658, f_geo=0.0387, f_vol=0.0639)
        return image.select([]).addBands([blue, green, red, nir, re1, re2, re3, nir, re4, swir1, swir2])

    def _correct_band(image, band_name, kvol, kvol0, f_iso, f_geo, f_vol):
        iso   = ee.Image(f_iso)
        geo   = ee.Image(f_geo)
        vol   = ee.Image(f_vol)
        pred  = vol.multiply(kvol).add(geo.multiply(kvol)).add(iso).rename(['pred'])
        pred0 = vol.multiply(kvol0).add(geo.multiply(kvol0)).add(iso).rename(['pred0'])
        cfac  = pred0.divide(pred).rename(['cfac'])
        corr  = image.select(band_name).multiply(cfac).rename([band_name])
        return corr

    def _kvol(sunAz, sunZen, viewAz, viewZen):
        relative_azimuth = sunAz.subtract(viewAz).rename(['relAz'])
        pa1 = viewZen.cos().multiply(sunZen.cos())
        pa2 = viewZen.sin().multiply(sunZen.sin()).multiply(relative_azimuth.cos())
        phase_angle1 = pa1.add(pa2)
        phase_angle  = phase_angle1.acos()
        p1 = ee.Image(PI.divide(2)).subtract(phase_angle)
        p2 = p1.multiply(phase_angle1)
        p3 = p2.add(phase_angle.sin())
        p4 = sunZen.cos().add(viewZen.cos())
        p5 = ee.Image(PI.divide(4))
        kvol = p3.divide(p4).subtract(p5).rename(['kvol'])
        viewZen0     = ee.Image(0)
        pa10         = viewZen0.cos().multiply(sunZen.cos())
        pa20         = viewZen0.sin().multiply(sunZen.sin()).multiply(relative_azimuth.cos())
        phase_angle10 = pa10.add(pa20)
        phase_angle0  = phase_angle10.acos()
        p10 = ee.Image(PI.divide(2)).subtract(phase_angle0)
        p20 = p10.multiply(phase_angle10)
        p30 = p20.add(phase_angle0.sin())
        p40 = sunZen.cos().add(viewZen0.cos())
        p50 = ee.Image(PI.divide(4))
        kvol0 = p30.divide(p40).subtract(p50).rename(['kvol0'])
        return [kvol, kvol0]

    def applyBRDF(image):
        date      = image.date()
        footprint = ee.List(image.geometry().bounds(3).bounds().coordinates().get(0))
        angles    = getsunAngles(date, footprint)
        sunAz     = angles[0]
        sunZen    = angles[1]
        viewAz  = azimuth(footprint)
        viewZen = zenith(footprint)
        kval  = _kvol(sunAz, sunZen, viewAz, viewZen)
        kvol  = kval[0]
        kvol0 = kval[1]
        result = _apply(image, kvol.multiply(PI), kvol0.multiply(PI))
        return result

    correctedImg = applyBRDF(img)
    return correctedImg


# =========================
# 7. COMPOSITING ==========
# =========================

# FUNCTION: quality_S2_CloudlessComposite
# USE: Produce cloud-free Sentinel-2 composites and export
# PARAMETERS:
#  dataset_id = numerical dataset identifier
#  feat_col = feature collection of buffered plot/site points corresponding to dataset_id
#  start = start date yyyy-mm-dd
#  end = end date yyyy-mm-dd
#  location = export location (asset vs. drive)
#  folder = if exporting to asset, specify folder
# NOTES: Parameters are populated automatically within 1_proc_site_representativeness script
# AUTHOR: Melissa Rose, NAU

def quality_S2_CloudlessComposite(dataset_id, feat_col, start, end, location, folder):

    # ----- INPUTS -----

    export_loc   = location
    asset_folder = folder
    aoi_name     = dataset_id
    aoi          = ee.Geometry(feat_col.geometry().bounds())
    start_date   = start
    end_date     = end

    # Sentinel-2 band selection
    save_bands = ['Green', 'Red', 'RE2', 'NIR', 'SWIR1', 'SWIR2',
                  'NDVI', 'EVI', 'RGVI', 'NDMI', 'NBR', 'SVVI']

    # Cloud filtering thresholds
    cloud_perc_filter = 70
    cloud_prob_thresh = 15
    nir_dark_thresh   = 0.15
    cloud_proj_dist   = 1
    buffer            = 50

    compos_meth = 'quality'

    # ----- PROCESSING ----

    s2_l2a_filt = (ee.ImageCollection('COPERNICUS/S2_SR')
                   .filterBounds(aoi)
                   .filter(ee.Filter.calendarRange(182, 228, 'day_of_year'))
                   .filterDate(start_date, end_date)
                   .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', cloud_perc_filter)))

    s2_cs_mask = L2A_mask_cloudless_wshadow(
        s2_l2a_filt, cloud_prob_thresh, nir_dark_thresh, cloud_proj_dist, buffer
    )

    s2_ind = (s2_cs_mask
              .map(scale_refl_toFlt)
              .map(sel_rename_bands)
              .map(brdf_corr)
              .map(add_spec_indices))

    # ----- COMPOSITING ----

    compos_index = 'NDVI'
    qual_perc    = 98

    s2_cm_index_perc = (ee.Image(s2_ind.select(compos_index)
                                 .reduce(ee.Reducer.percentile([qual_perc])))
                        .rename('ref_val'))

    def diff_index_perc(img):
        diff = (ee.Image.constant(100)
                .subtract(img.select([compos_index])
                            .subtract(s2_cm_index_perc.select(['ref_val']))
                            .abs())
                .rename('qual_val'))
        return img.addBands(diff)

    s2_qm = (s2_ind
             .map(diff_index_perc)
             .qualityMosaic('qual_val')
             .select(save_bands)
             .clip(aoi))

    s2_compos = s2_qm.multiply(10000).int16()

    # ----- OUTPUTS -----

    # Export options are commented out to preserve the original JS logic.
    # Uncomment and adjust as needed for Python batch export.

    # asset_name = f'{aoi_name}_S2_L2A_compos_{compos_meth}_{start_date}_{end_date}'
    # if export_loc == 'asset':
    #     ee.batch.Export.image.toAsset(
    #         image=s2_compos, assetId=asset_folder + asset_name,
    #         description=asset_name, region=aoi, scale=10,
    #         crs='EPSG:4326', maxPixels=1e13
    #     ).start()
    # elif export_loc == 'drive':
    #     ee.batch.Export.image.toDrive(
    #         image=s2_compos, description=asset_name,
    #         fileNamePrefix=asset_name, region=aoi, scale=10,
    #         crs='EPSG:4326', maxPixels=1e13,
    #         skipEmptyTiles=True, fileFormat='GeoTIFF'
    #     ).start()

    return s2_compos
