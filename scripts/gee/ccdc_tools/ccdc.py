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

Utility functions for filtering and modifying CCD output.
"""

import math
import ee


# =========================
# 1. FUNCTIONS ==============
# =========================

def buildSegmentTag(n_segments):
    """
    Create sequence of segment strings.
    @param {Integer} n_segments  Number of segments to create labels for
    @returns {ee.List} List of segment names (e.g. S1, S2)
    """
    return ee.List.sequence(1, n_segments).map(lambda i: ee.String('S').cat(ee.Number(i).int()))


def buildBandTag(tag, band_list):
    """
    Create sequence of band names for a given string tag.
    @param {string} tag       String tag to use (e.g. 'RMSE')
    @param {array}  band_list List of band names to combine with tag
    @returns {ee.List} List of band names combined with tag name
    """
    bands = ee.List(band_list)
    return bands.map(lambda s: ee.String(s).cat('_' + tag))


def buildMagnitude(fit, n_segments, band_list):
    """
    Extract CCDC magnitude image from current CCDC result format.
    @param {ee.Image} fit         Image with CCD results
    @param {number}   n_segments  Number of segments to extract
    @param {array}    band_list   Client-side list with band names to use
    @returns {ee.Image} Image with magnitude of change per segment per band
    """
    segment_tag = buildSegmentTag(n_segments)
    zeros = ee.Image(ee.Array(ee.List.repeat(0, n_segments)))

    def retrieve_mags(band):
        # Pad zeroes for pixels that have less than n_segments, then slice the first n_segments values
        mag_img = fit.select(band + '_magnitude').arrayCat(zeros, 0).float().arraySlice(0, 0, n_segments)
        tags = segment_tag.map(lambda x: ee.String(x).cat('_').cat(band).cat('_MAG'))
        return mag_img.arrayFlatten([tags])

    return ee.Image(band_list.map(retrieve_mags))


def buildRMSE(fit, n_segments, band_list):
    """
    Extract CCDC RMSE image from current CCDC formatted results.
    @param {ee.Image} fit         Image with CCDC results
    @param {number}   n_segments  Number of segments to extract
    @param {array}    band_list   Client-side list with band names to use
    @returns {ee.Image} Image with RMSE of each segment per band
    """
    segment_tag = buildSegmentTag(n_segments)
    zeros = ee.Image(ee.Array(ee.List.repeat(0, n_segments)))

    def retrieve_mags(band):
        # Pad zeroes for pixels that have less than n_segments, then slice the first n_segments values
        mag_img = fit.select(band + '_rmse').arrayCat(zeros, 0).float().arraySlice(0, 0, n_segments)
        tags = segment_tag.map(lambda x: ee.String(x).cat('_').cat(band).cat('_RMSE'))
        return mag_img.arrayFlatten([tags])

    return ee.Image(band_list.map(retrieve_mags))


def buildCoefs(fit, n_segments, band_list):
    """
    Extract CCDC Coefficients from current CCDC formatted result.
    @param {ee.Image} fit         Image with CCD results
    @param {number}   n_segments  Number of segments to extract
    @param {array}    band_list   Client-side list with band names to use
    @returns {ee.Image} Image with coefficients per band
    """
    segment_tag = buildSegmentTag(n_segments)
    harmonic_tag = ['INTP', 'SLP', 'COS', 'SIN', 'COS2', 'SIN2', 'COS3', 'SIN3']
    zeros = ee.Image(ee.Array([ee.List.repeat(0, len(harmonic_tag))])).arrayRepeat(0, n_segments)

    def retrieve_coefs(band):
        coef_img = fit.select(band + '_coefs').arrayCat(zeros, 0).float().arraySlice(0, 0, n_segments)
        tags = segment_tag.map(lambda x: ee.String(x).cat('_').cat(band).cat('_coef'))
        return coef_img.arrayFlatten([tags, harmonic_tag])

    return ee.Image(band_list.map(retrieve_coefs))


def buildStartEndBreakProb(fit, n_segments, tag):
    """
    Extract data for CCDC 1D-array, non-spectral bands
    (tStart, tEnd, tBreak, changeProb or numObs).
    @param {ee.Image} fit         Image with CCD results
    @param {integer}  n_segments  Number of segments to extract
    @param {string}   tag         Client-side string to use as name in the output bands
    @returns {ee.Image} Image with values for tStart, tEnd, tBreak, changeProb or numObs
    """
    segment_tag = buildSegmentTag(n_segments).map(lambda s: ee.String(s).cat('_' + tag))
    zeros = ee.Array(0).repeat(0, n_segments)
    mag_img = fit.select(tag).arrayCat(zeros, 0).float().arraySlice(0, 0, n_segments)
    return mag_img.arrayFlatten([segment_tag])


def buildCcdImage(fit, n_segments, band_list):
    """
    Transform ccd results from array image to 'long' multiband format.
    @param {ee.Image} fit         Image with CCD results
    @param {number}   n_segments  Number of segments to extract
    @param {array}    band_list   Client-side list with band names to use
    @returns {ee.Image} Image with all results from CCD in 'long' image format
    """
    magnitude = buildMagnitude(fit, n_segments, band_list)
    rmse = buildRMSE(fit, n_segments, band_list)
    coef = buildCoefs(fit, n_segments, band_list)
    t_start = buildStartEndBreakProb(fit, n_segments, 'tStart')
    t_end = buildStartEndBreakProb(fit, n_segments, 'tEnd')
    t_break = buildStartEndBreakProb(fit, n_segments, 'tBreak')
    probs = buildStartEndBreakProb(fit, n_segments, 'changeProb')
    nobs = buildStartEndBreakProb(fit, n_segments, 'numObs')
    return ee.Image.cat(coef, rmse, magnitude, t_start, t_end, t_break, probs, nobs)


def getSyntheticForYear(image, date, date_format, band, segs):
    """
    Create synthetic image for specified band.
    @param {ee.Image} image        Image with CCD results in long multi-band format
    @param {number}   date         Date to extract the segments for, in the format that ccd was run in
    @param {number}   date_format  Code of the date format that ccdc was run in (e.g. 1 for frac years)
    @param {string}   band         Band name to use for creation of synthetic image
    @param {array}    segs         List of segment names to use
    @returns {ee.Image} Synthetic image for the given date and band
    """
    tfit = date
    PI2 = 2.0 * math.pi
    OMEGAS = [PI2 / 365.25, PI2, PI2 / (1000 * 60 * 60 * 24 * 365.25)]
    omega = OMEGAS[date_format]
    image_t = ee.Image.constant([
        1, tfit,
        tfit.multiply(omega).cos(),
        tfit.multiply(omega).sin(),
        tfit.multiply(omega * 2).cos(),
        tfit.multiply(omega * 2).sin(),
        tfit.multiply(omega * 3).cos(),
        tfit.multiply(omega * 3).sin()
    ]).float()

    COEFS = ['INTP', 'SLP', 'COS', 'SIN', 'COS2', 'SIN2', 'COS3', 'SIN3']
    new_params = getMultiCoefs(image, date, [band], COEFS, False, segs, 'before')
    return image_t.multiply(new_params).reduce('sum').rename(band)


def getMultiSynthetic(image, date, date_format, band_list, segs):
    """
    Create synthetic image for a list of bands.
    @param {ee.Image} image        Image with CCD results in long multi-band format
    @param {number}   date         Date to extract the segments for, in the format that ccd was run in
    @param {number}   date_format  Code of the date format that ccdc was run in (e.g. 1 for frac years)
    @param {array}    band_list    List of bands to get synthetic data for
    @param {array}    segs         List of segment names to use
    @returns {ee.Image} Synthetic image for the given date and bands
    """
    def retrieve_synthetic(band):
        return getSyntheticForYear(image, date, date_format, band, segs)
    return ee.Image.cat(band_list.map(retrieve_synthetic))


def fillNoData(fit, n_coefs, n_bands):
    """
    @deprecated (No longer necessary)
    Replace nodata in CCD output and fill with zeros.
    @param {ee.Image} fit      Image with CCD results
    @param {number}   n_coefs  Number of coefficients present in the results
    @param {number}   n_bands  Number of spectral bands used to produce the results
    @returns {ee.Image} Image with nodata areas replaced with zeros
    """
    d1 = ee.Image(ee.Array([0]).double())
    d2 = ee.Image(ee.Array([ee.List.repeat(-9999, n_coefs)])).double()
    upper = ee.Image([d1, d1, d1, d1.int32(), d1])

    arr_center = []
    arr_bottom = []
    for i in range(n_bands):
        arr_center.append(d2)
        arr_bottom.extend([d1, d1])
    center = ee.Image(arr_center)
    bottom = ee.Image(arr_bottom)

    mock = (upper.addBands(center).addBands(bottom)
                 .rename(fit.bandNames()).updateMask(fit.mask()))
    new_image = ee.ImageCollection([mock, fit]).mosaic()
    return new_image


def dateToDays(str_date):
    """
    @deprecated (No longer necessary)
    Return a date as days from 01-01-0000.
    @param {String} str_date  Date in the format accepted by ee.Date
    @returns {ee.Number} Date expressed as days since 01-01-0000
    """
    date = ee.Date(str_date)
    # Number of days since 01-01-0000 until 01-01-1970
    epoch = ee.Number(719177)
    # Convert millis to days
    days = ee.Number(date.millis().divide(86400000))
    return days.add(epoch)


def dateToSegment(ccd_results, date, seg_names):
    """
    @deprecated (No longer necessary) Find segments that intersect a given date.
    @param {ee.Image}  ccd_results  CCD results in long multi-band format
    @param {Number}    date         Date in the format that was used to run CCD
    @param {ee.List}   seg_names    Segment names to use
    @returns {ee.Image} segmentMatch  Mask image indicating segments that intersect the given date
    """
    start_bands = ccd_results.select('.*_tStart').rename(seg_names)
    end_bands = ccd_results.select('.*_tEnd').rename(seg_names)
    start = start_bands.lte(date)
    end = end_bands.gte(date)
    segment_match = start.And(end)
    return segment_match


def filterCoefs(ccd_results, date, band, coef, seg_names, behavior):
    """
    Filter coefficients for a given date using a mask.
    @param {ee.Image} ccd_results  CCD results in long multi-band format
    @param {string}   date         Date in the same format that CCD was run with
    @param {string}   band         Band to select
    @param {string}   coef         Coef to select. Options are "INTP", "SLP", "COS", "SIN",
                                   "COS2", "SIN2", "COS3", "SIN3", "RMSE", "MAG"
    @param {ee.List}  seg_names    List of segment names to use
    @param {String}   behavior     Method to find intersecting ('normal') or closest
                                   segment to given date ('before' or 'after')
    @returns {ee.Image} Single band image with the values for the selected band/coefficient
    """
    start_bands = ccd_results.select('.*_tStart').rename(seg_names)
    end_bands = ccd_results.select('.*_tEnd').rename(seg_names)

    # Get all segments for a given band/coef. Underscore in concat ensures
    # that bands with similar names are not selected twice (e.g. GREEN, GREENNESS)
    sel_str = '.*' + band + '_.*' + coef  # Client side concat
    coef_bands = ccd_results.select(sel_str)

    # Select a segment based on conditions
    if behavior == 'normal':
        start = start_bands.lte(date)
        end = end_bands.gte(date)
        segment_match = start.And(end)
        out_coef = coef_bands.updateMask(segment_match).reduce(ee.Reducer.firstNonNull())
    elif behavior == 'after':
        segment_match = end_bands.gt(date)
        out_coef = coef_bands.updateMask(segment_match).reduce(ee.Reducer.firstNonNull())
    elif behavior == 'before':
        # Mask start to avoid comparing against zero, mask after to remove zeros from logical comparison
        segment_match = start_bands.selfMask().lt(date).selfMask()
        out_coef = coef_bands.updateMask(segment_match).reduce(ee.Reducer.lastNonNull())
    else:
        out_coef = coef_bands.reduce(ee.Reducer.firstNonNull())

    # TODO: Add an "automatic" after, then before behavior
    return out_coef


def normalizeIntercept(intercept, start, end, slope):
    """
    Normalize the intercept to the middle of the segment time period,
    instead of the 0 time period.
    @param {ee.Image} intercept  Image band representing model intercept
    @param {ee.Image} start      Image band representing model slope date
    @param {ee.Image} end        Image band representing model end date
    @param {ee.Image} slope      Image band representing model slope
    @returns {ee.Image} Image band representing normalized intercept
    """
    middle_date = ee.Image(start).add(ee.Image(end)).divide(2)
    slope_coef = ee.Image(slope).multiply(middle_date)
    return ee.Image(intercept).add(slope_coef)


def getCoef(ccd_results, date, band_list, coef, seg_names, behavior):
    """
    Get image with a single coefficient for all bands.
    @param {ee.Image} ccd_results  CCD results in long multi-band format
    @param {string}   date         Date in the same format that CCD was run with
    @param {array}    band_list    List of all bands to include
    @param {array}    coef         Coef to select. Options are "INTP", "SLP", "COS", "SIN",
                                   "COS2", "SIN2", "COS3", "SIN3", "RMSE", "MAG"
    @param {ee.List}  seg_names    List of segment names to use
    @param {string}   behavior     Method to find intersecting ('normal') or closest
                                   segment to given date ('before' or 'after')
    @returns {ee.Image} coefs  Image with the values for the selected bands x coefficient
    """
    def inner(band):
        band_coef = filterCoefs(ccd_results, date, band, coef, seg_names, behavior)
        return band_coef.rename(band + '_' + coef)  # Client side concat
    coefs = ee.Image(band_list.map(inner))  # Client side map
    return coefs


def applyNorm(band_coefs, seg_start, seg_end):
    """
    Apply normalization to intercepts.
    @param {ee.Image} band_coefs  Band x coefficients image. Must include slopes
    @param {ee.Image} seg_start   Image with dates representing the start of the segment
    @param {ee.Image} seg_end     Image with dates representing the end of the segment
    @returns {ee.Image} band_coefs  Updated input image with normalized intercepts
    """
    intercepts = band_coefs.select('.*INTP')
    slopes = band_coefs.select('.*SLP')
    normalized = normalizeIntercept(intercepts, seg_start, seg_end, slopes)
    return band_coefs.addBands(**{'srcImg': normalized, 'overwrite': True})


def getMultiCoefs(ccd_results, date, band_list, coef_list, cond, seg_names, behavior):
    """
    Get image with bands x coefficients given in a list.
    @param {ee.Image} ccd_results  CCD results in long multi-band format
    @param {string}   date         Date in the same format that CCD was run with
    @param {array}    band_list    List of all bands to include
    @param {list}     coef_list    List of coefs to select. Options are "INTP", "SLP",
                                   "COS", "SIN", "COS2", "SIN2", "COS3", "SIN3", "RMSE", "MAG"
    @param {boolean}  cond         Normalize intercepts? If True, requires "INTP" and "SLP"
                                   to be selected in coef_list
    @param {ee.List}  seg_names    List of segment names to use
    @param {string}   behavior     Method to find intersecting ('normal') or closest
                                   segment to given date ('before' or 'after')
    @returns {ee.Image} coefs  Image with the values for the selected bands x coefficients
    """
    def inner(coef):
        return getCoef(ccd_results, date, band_list, coef, seg_names, behavior)

    coefs = ee.Image(coef_list.map(inner))

    # Normalized
    seg_start = filterCoefs(ccd_results, date, '', 'tStart', seg_names, behavior)
    seg_end = filterCoefs(ccd_results, date, '', 'tEnd', seg_names, behavior)
    norm_coefs = applyNorm(coefs, seg_start, seg_end)

    out_coefs = ee.Algorithms.If(cond, norm_coefs, coefs)
    return ee.Image(out_coefs)


def getChanges(ccd_results, start_date, end_date, seg_names):
    """
    Filter segments with change in a given range.
    @param {ee.Image} ccd_results  CCD results in long multi-band format
    @param {Number}   start_date   Start date in the format that was used to run CCD
    @param {Number}   end_date     End date in the format that was used to run CCD
    @param {ee.List}  seg_names    List of segment names matching the number of segments in the bands
    @returns {ee.Image} Mask image indicating which pixel/segments have changes in the specified time range
    """
    break_bands = ccd_results.select('.*_tBreak').rename(seg_names)
    segment_match = break_bands.gte(start_date).And(break_bands.lt(end_date))
    return segment_match


def filterMag(ccd_results, start_date, end_date, band, seg_names):
    """
    Obtain change with largest magnitude, timing of that break, and total number of breaks
    for a given date range and band.
    @param {ee.Image} ccd_results  CCD results in long multi-band format
    @param {number}   start_date   Start date in the format that was used to run CCD
    @param {number}   end_date     End date in the format that was used to run CCD
    @param {string}   band         Spectral band
    @param {ee.List}  seg_names    List of segment names matching the number of segments in the bands
    @returns {ee.Image} Image with three bands: magnitude of largest break, timing, and total number of breaks
    """
    seg_mask = getChanges(ccd_results, start_date, end_date, seg_names)
    sel_str = '.*' + band + '.*' + 'MAG'  # Client side concat
    feat_bands = ccd_results.select(sel_str)

    # Need abs vals because mags can be negative too!
    filtered_mag = feat_bands.mask(seg_mask)
    filtered_abs_mag = filtered_mag.abs()
    max_abs_mag = filtered_abs_mag.reduce(ee.Reducer.max())

    # Find which 'index' matches that abs mag
    matched_mag_mask = filtered_abs_mag.eq(max_abs_mag)

    # Use that index to select the magnitude with the original sign, and the timing of that break
    selected_mag = filtered_mag.mask(matched_mag_mask).reduce(ee.Reducer.firstNonNull())
    filtered_t_break = ccd_results.select('.*tBreak').mask(matched_mag_mask).reduce(ee.Reducer.firstNonNull())
    num_t_break = ccd_results.select('.*tBreak').mask(seg_mask).reduce(ee.Reducer.count())

    return (selected_mag.addBands(filtered_t_break)
                        .addBands(num_t_break)
                        .rename(['MAG', 'tBreak', 'numTbreak']))


def phaseAmplitude(img, bands, sin_name, cos_name):
    """
    Get phase and amplitude for a single spectral band.
    @param {ee.Image} img        CCD results in long multi-band format
    @param {List}     bands      List with the name of the bands for which to calculate ampl. and phase
    @param {String}   sin_name   Band suffix of the desired sine harmonic coefficient (e.g. '_SIN')
    @param {String}   cos_name   Band suffix of the desired cosine harmonic coefficient (e.g. '_COS')
    @returns {ee.Image} Image with two bands representing phase and amplitude of the desired harmonic
    """
    sin_names = bands.map(lambda x: x + sin_name)
    cos_names = bands.map(lambda x: x + cos_name)
    phase_names = bands.map(lambda x: x + '_PHASE')
    amplitude_names = bands.map(lambda x: x + '_AMPLITUDE')
    phase = (img.select(sin_names).atan2(img.select(cos_names))
                # Scale to [0, 1] from radians.
                .unitScale(-math.pi, math.pi)
                .multiply(365)  # To get phase in days!
                .rename(phase_names))
    amplitude = img.select(sin_names).hypot(img.select(cos_names)).rename(amplitude_names)
    return phase.addBands(amplitude)


def newPhaseAmplitude(img, sin_expr, cos_expr):
    """
    Get phase and amplitude. Replace old function with this.
    @param {ee.Image} img       CCD results in long multi-band format
    @param {String}   sin_expr  Regular expression of the sine harmonic coefficient
                                (e.g. '.*SIN.*') for all harmonics
    @param {String}   cos_expr  Regular expression of the cosine harmonic coefficient
                                (e.g. '.*COS.*') for all harmonics.
                                Must retrieve the same number of bands as sin_expr
    @returns {ee.Image} Image with bands representing phase and amplitude of the desired harmonic
    """
    sin = img.select(sin_expr)
    cos = img.select(cos_expr)

    phase = (sin.atan2(cos)
               # Scale to [0, 1] from radians.
               .unitScale(-3.14159265359, 3.14159265359)
               .multiply(365))  # To get phase in days!

    amplitude = sin.hypot(cos)

    phase_names = phase.bandNames().map(lambda x: ee.String(x).replace('_SIN', '_PHASE'))
    amplitude_names = amplitude.bandNames().map(lambda x: ee.String(x).replace('_SIN', '_AMPLITUDE'))
    return phase.rename(phase_names).addBands(amplitude.rename(amplitude_names))
