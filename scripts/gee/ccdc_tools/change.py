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

Utility functions for exploring and processing change information.
"""

import math
import ee
from . import dates as date_utils
from . import stats as stats_utils


# =========================
# 1. FUNCTIONS ==============
# =========================

def createTimeBand(img):
    """
    Add time band and constant to an image.
    @param {ee.Image} img  Input image
    @returns {ee.Image} Image with added 't' and 'constant' bands
    """
    year = ee.Image(img.date().difference('1970-01-01', 'year')).rename('t')
    constant = ee.Image(1).rename('constant')
    return img.addBands(constant).addBands(year.float())


def constructBandNames(base, lst):
    """
    Get a sequence of names using a base name and a list of strings to append.
    @param {string}   base  Base name string
    @param {ee.List}  lst   List of values to append
    @returns {ee.List} List of concatenated band names
    """
    return ee.List(lst).map(lambda i: ee.String(base).cat(ee.Number(i).int()))


def addDependents(image):
    """
    Add a constant and time band to an image.
    @param {ee.Image} image  Input image
    @returns {ee.Image} Image with added 'constant' and 't' bands
    """
    # Compute time in fractional years since the epoch.
    years = image.date().difference('1970-01-01', 'year')
    time_radians = ee.Image(years.multiply(2 * math.pi)).rename('t')
    constant = ee.Image(1)
    return image.addBands(constant).addBands(time_radians.float())


def calc_scaled_cusum(resids, crit_val):
    """
    Calculate scaled CUSUM OLS residuals.
    @param {ee.ImageCollection} resids    Collection of residual images
    @param {number}             crit_val  Critical value threshold for break detection
    @returns {ee.ImageCollection} Scaled CUSUM collection with 'cusum' and 'break_detected' bands
    """
    # Initialize empty cusum
    dt = ee.Date(resids.first().get('system:time_start')).advance(-1, 'day').millis()
    cs_0 = ee.List([ee.Image(0.0001).set('system:time_start', dt).rename('cusum')])

    # Simple CUSUM function. Unmask images to prevent sum to fail due to masked values.
    # Unmasked values will have a value of zero and won't affect the calculation
    def cusum_calc(image, lst):
        previous = ee.Image(ee.List(lst).get(-1))
        current_cusum = (previous.select('cusum').unmask()
                                 .add(image.unmask())
                                 .set('system:time_start', image.get('system:time_start')))
        return ee.List(lst).add(current_cusum.rename('cusum'))

    # Calculate CUSUM
    cusum = ee.ImageCollection(ee.List(resids.iterate(cusum_calc, cs_0)))

    # Scaled cumulative residuals, based on statsmodels (Ploberger, Werner, and Walter Kramer 1992)
    nobs = resids.count()  # Count excludes masked vals, which is what we want here
    nobs_sigma2 = resids.map(lambda x: x.pow(2)).sum()
    ddof = ee.Number(8)  # # of Harmonic regression pairs + intp + slope. Automate this calculation
    # This notation replicates the python operator order
    nobs_sigma2 = nobs_sigma2.divide(nobs.subtract(ddof)).multiply(nobs)

    # Actual calculation of scaled CUSUM. Append breakpoint band
    def scale_image(img):
        scaled_res = img.divide(nobs_sigma2.sqrt()).set('system:time_start', img.get('system:time_start'))
        breakpoint = scaled_res.abs().gte(crit_val).rename('break_detected')
        return scaled_res.addBands(breakpoint)

    scaled_resid = cusum.map(scale_image)
    return scaled_resid


def addFracyearBand(img):
    """
    Add date band in fractional year format to compare against tstart/tend.
    @param {ee.Image} img  Input image
    @returns {ee.Image} Image with added 'fracYear' band
    """
    return img.addBands(ee.Image(date_utils.msToFrac(img.date().millis())).rename('fracYear'))


def harmonic_regression(ts_collection, band):
    """
    Run a harmonic regression for a given band in a collection and return
    observed and predicted values, residuals, and coefficients.
    @param {ee.ImageCollection} ts_collection  Time series image collection
    @param {string}             band           Band name to regress
    @returns {ee.List} [harmonicTrendCoefficients, regression_results]
    """
    # The number of cycles per year to model and variable name to use
    harmonics = 3
    dependent = band

    # Make a list of harmonic frequencies to model.
    harmonic_frequencies = ee.List.sequence(1, harmonics)

    # Construct lists of names for the harmonic terms.
    cos_names = constructBandNames('cos_', harmonic_frequencies)
    sin_names = constructBandNames('sin_', harmonic_frequencies)

    # Independent variables.
    independents = ee.List(['constant', 't']).cat(cos_names).cat(sin_names)

    # Function to compute the specified number of harmonics and add them as bands.
    # Assumes the time band is present.
    def add_harmonics(freqs):
        def inner(image):
            # Make an image of frequencies.
            frequencies = ee.Image.constant(freqs)
            # This band should represent time in radians.
            time = ee.Image(image).select('t')
            # Get the cosine terms.
            cosines = time.multiply(frequencies).cos().rename(cos_names)
            # Get the sin terms.
            sines = time.multiply(frequencies).sin().rename(sin_names)
            return image.addBands(cosines).addBands(sines)
        return inner

    # Add variables.
    harmonic_landsat = ts_collection.map(addDependents).map(add_harmonics(harmonic_frequencies))

    # The output of the regression reduction is a 4x1 array image.
    harmonic_trend = (harmonic_landsat
                      .select(independents.add(dependent))
                      .reduce(ee.Reducer.linearRegression(independents.length(), 1)))

    # Turn the array image into a multi-band image of coefficients.
    harmonic_trend_coefficients = (harmonic_trend.select('coefficients')
                                                  .arrayProject([0])
                                                  .arrayFlatten([independents]))

    # Compute fitted values.
    def compute_fitted(image):
        return image.addBands(
            image.select(independents)
                 .multiply(harmonic_trend_coefficients)
                 .reduce('sum')
                 .rename('fitted'))
    fitted_harmonic = harmonic_landsat.map(compute_fitted)

    # Compute residuals manually
    def calc_residuals(img):
        resid = img.select(dependent).subtract(img.select('fitted')).rename('residuals')
        return img.addBands(resid)
    regression_results = fitted_harmonic.map(calc_residuals)

    return ee.List([harmonic_trend_coefficients, regression_results])


def getOmission(landsat_col, ccdc_img, seg_id, band, crit_val):
    """
    Test integration with CCDC segments — omission test.
    @param {ee.ImageCollection} landsat_col  Cloud and shadow filtered Landsat collection
    @param {ee.Image}           ccdc_img     CCDC results in long format
    @param {string}             seg_id       Segment ID to apply function to (e.g. 'S1', 'S2')
    @param {string}             band         Spectral band to calculate test on (e.g. 'SWIR1')
    @param {number}             crit_val     Critical value for CUSUM threshold
    @returns {ee.Image} Image with bands 'omission' and 'maxCusum'
    """
    # a) Get tstart and tend for current segment
    t_start = ccdc_img.select(seg_id + '_tStart')
    t_end = ccdc_img.select(seg_id + '_tEnd')

    # Add band with fractional year
    frac_year_col = landsat_col.map(addFracyearBand)

    # b) Compare date for each image against the tstart-tend per pixel, and mask accordingly
    def mask_dates(img):
        frac_year = img.select('fracYear')
        mask = ee.Image(frac_year.gte(t_start)).And(ee.Image(frac_year.lte(t_end))).rename('dateMask')
        return img.updateMask(mask)

    masked_col = frac_year_col.map(mask_dates)
    # Do harmonic regression instead of simple linear
    regression_results = ee.ImageCollection(harmonic_regression(masked_col, band).get(1))

    # f) Calculate scaled residuals and occurrence of break if threshold reached
    residuals = regression_results.select('residuals')
    # Calculate scaled CUSUM and time of break for each image
    scaled_cusum = calc_scaled_cusum(residuals, crit_val)
    cusum = scaled_cusum.select('cusum').map(lambda x: x.abs())
    break_detected = scaled_cusum.select('break_detected')
    # Get presence of any breaks detected for the current segment, and the max abs value
    cusum_break = break_detected.reduce(ee.Reducer.anyNonZero()).rename('omission')
    cusum_threshold = cusum.reduce(ee.Reducer.max()).rename('maxCusum')

    return ee.Image.cat([cusum_break, cusum_threshold])


def getCommission(landsat_col, ccdc_img, seg1, seg2, band):
    """
    Test integration with CCDC segments — commission test (Chow test).
    @param {ee.ImageCollection} landsat_col  Cloud and shadow filtered Landsat collection
    @param {ee.Image}           ccdc_img     CCDC results in long format
    @param {string}             seg1         First segment ID (e.g. 'S1')
    @param {string}             seg2         Second segment ID (e.g. 'S2')
    @param {string}             band         Spectral band to calculate test on (e.g. 'SWIR1')
    @returns {ee.Image} Image with bands 'chow_F' and 'prob'
    """
    # a) Get tstart and tend for each segment
    t_start1 = ccdc_img.select(seg1 + '_tStart')
    t_end1 = ccdc_img.select(seg1 + '_tEnd')
    t_start2 = ccdc_img.select(seg2 + '_tStart')
    t_end2 = ccdc_img.select(seg2 + '_tEnd')

    # Add band with fractional year
    frac_year_col = landsat_col.map(addFracyearBand)

    # b) Compare date for each image against the tstart-tend per pixel, and mask accordingly
    def mask_dates(col, t_start, t_end):
        def inner(img):
            frac_year = img.select('fracYear')
            mask = ee.Image(frac_year.gte(t_start)).And(ee.Image(frac_year.lte(t_end))).rename('dateMask')
            return img.updateMask(mask)
        return col.map(inner)

    masked_col1 = mask_dates(frac_year_col, t_start1, t_end1)
    masked_col2 = mask_dates(frac_year_col, t_start2, t_end2)
    masked_col_global = mask_dates(frac_year_col, t_start1, t_end2)

    # Run harmonic regression instead of linear
    regression_results1 = ee.ImageCollection(harmonic_regression(masked_col1, band).get(1))
    regression_results2 = ee.ImageCollection(harmonic_regression(masked_col2, band).get(1))
    regression_results_global = ee.ImageCollection(harmonic_regression(masked_col_global, band).get(1))

    resid1 = regression_results1.select('residuals')
    resid2 = regression_results2.select('residuals')
    resid_global = regression_results_global.select('residuals')

    # 6) Get sum of squared residuals for each model
    def sum_squares(col):
        return col.map(lambda x: x.pow(2)).sum()

    ss1 = sum_squares(resid1)
    ss2 = sum_squares(resid2)
    ss_gl = sum_squares(resid_global)

    # Get nobs per segment, directly from CCDC record
    n1 = ccdc_img.select(seg1 + '_numObs')
    n2 = ccdc_img.select(seg2 + '_numObs')
    n_gl = n1.add(n2)

    # 7) Calculate Chow test statistic
    k = ee.Number(4)  # Num of parameters, e.g. 2 pairs of harmonics
    num = ss_gl.subtract(ss1.add(ss2)).divide(k.add(1))
    den = ss1.add(ss2).divide(n_gl.subtract(k.add(1).multiply(2)))
    chow = num.divide(den)

    # Calculate probability
    f_prob = stats_utils.F_cdf(chow, k, n_gl.subtract(k.multiply(2)))
    return chow.addBands(f_prob).rename(['chow_F', 'prob'])
