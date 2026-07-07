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

Date functions for converting from and to the different date formats
used by the CCDC algorithm.
"""

import ee

# =========================
# 1. CONSTANTS ==============
# =========================

# Conversion factor from ms to days
MS_TO_DAYS = 86400000
# Number of days in common era until epoch 01-01-1970 (non-inclusive)
EPOCH_DAYS = ee.Number(719529)


# =========================
# 2. FUNCTIONS ==============
# =========================

def msToDays(ms):
    """
    Convert milliseconds since epoch (01-01-1970) to number of days.
    @param {Number} ms  ms since 01-01-1970
    @returns {ee.Number} number of days
    """
    return ee.Number(ms).divide(MS_TO_DAYS)


def dateToJdays(str_date):
    """
    Convert Date to julian days in common era (i.e. days since 00-00-0000).
    @param {String} str_date  Date string in yyyy-mm-dd format
    @returns {ee.Number} Julian day
    """
    if not str_date:
        return 'Required parameter [str_date] missing'
    date = ee.Date(str_date)
    # Convert unix time to days
    return msToDays(date.millis()).add(EPOCH_DAYS)


def jdaysToms(jdays):
    """
    Convert julian day in common era (i.e. days since 00-00-0000) to ms since 1970-01-01.
    @param {Number} jdays  Julian day
    @returns {ee.Number} ms since 1970-01-01
    """
    day_since_epoch = ee.Number(jdays).subtract(EPOCH_DAYS)
    return day_since_epoch.multiply(MS_TO_DAYS)


def jdaysToDate(jdays):
    """
    Convert julian day in common era (i.e. days since 00-00-0000) to ee.Date.
    @param {Number} jdays  Julian day
    @returns {ee.Date} ee.Date
    """
    return ee.Date(jdaysToms(jdays))


def msToJdays(ms):
    """
    Convert ms since 1970-01-01 to julian day in common era (i.e. days since 00-00-0000).
    @param {Number} ms  ms since 1970-01-01
    @returns {ee.Number} Julian day
    """
    return ee.Number(msToDays(ms)).add(EPOCH_DAYS)


def msToFrac(ms):
    """
    Convert ms since 1970-01-01 to fractional year.
    @param {Number} ms  ms since 1970-01-01
    @returns {ee.Number} Fractional year
    """
    year = ee.Date(ms).get('year')
    frac = ee.Date(ms).getFraction('year')
    return year.add(frac)


def fracToms(frac):
    """
    Convert fractional time to ms since 1970-01-01.
    @param {Number} frac  Fractional year
    @returns {ee.Number} ms since 1970-01-01
    """
    fyear = ee.Number(frac)
    year = fyear.floor()
    d = fyear.subtract(year).multiply(365.25)
    day_one = ee.Date.fromYMD(year, 1, 1)
    return day_one.advance(d, 'day').millis()


def fracToDate(frac):
    """
    Convert fractional time to ee.Date.
    @param {Number} frac  Fractional year
    @returns {ee.Date} date object
    """
    ms = fracToms(frac)
    return msToDate(ms)


def msToDate(ms):
    """
    Convert ms to ee.Date.
    @param {number} ms  jdays as milliseconds
    @returns {ee.Date} ee.Date
    """
    return jdaysToDate(msToJdays(ms))


def convertDate(options):
    """
    Convert between any two date formats.

    @param {dict} options  parameter dictionary with keys:
        inputFormat  {Number}  date format according to ee ccdc function or 3 for string
                               0 = julian days, 1 = fractional year, 2 = ms, 3 = date string
        inputDate    {Object}  date as Number or String format matching inputFormat
        outputFormat {Number}  output date format according to ee ccdc function or 4 for ee.Date
                               0 = julian days, 1 = fractional year, 2 = ms, 4 = ee.Date
    @returns {Object} reformatted date
    """
    input_format = (options and options.get('inputFormat')) or 0
    input_date = (options and options.get('inputDate')) or None
    output_format = (options and options.get('outputFormat')) or 0

    if not input_date:
        return 'Required parameter [inputDate] missing'

    # First convert to millis
    if input_format == 0:
        milli = jdaysToms(input_date)
    elif input_format == 1:
        milli = fracToms(input_date)
    elif input_format == 2:
        milli = input_date
    elif input_format == 3:
        milli = jdaysToms(dateToJdays(input_date))
    else:
        milli = jdaysToms(input_date)

    # Now convert to output format
    if output_format == 0:
        output = msToJdays(milli)
    elif output_format == 1:
        output = msToFrac(milli)
    elif output_format == 2:
        output = milli
    elif output_format == 4:
        output = jdaysToDate(msToJdays(milli))
    else:
        output = msToJdays(milli)

    return output
