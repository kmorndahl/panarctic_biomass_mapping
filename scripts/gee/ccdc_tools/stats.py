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

Functions related to statistical distributions and tests.
"""

import ee

# =========================
# 1. CONSTANTS ==============
# =========================

UNIT = ee.Number(1)


# =========================
# 2. FUNCTIONS ==============
# =========================

def betaNum(x, y):
    """
    Beta function (scalar / ee.Number version).
    @param {Number} x  first parameter
    @param {Number} y  second parameter
    @returns {ee.Number} beta(x, y)
    """
    x = ee.Number(x)
    y = ee.Number(y)
    num = x.gamma().multiply(y.gamma())
    den = x.add(y).gamma()
    return ee.Number(num.divide(den))


def beta(x, y):
    """
    Beta function (image version).
    @param {ee.Image} x  first parameter image
    @param {ee.Image} y  second parameter image
    @returns {ee.Image} beta(x, y) image
    """
    x = ee.Image(x)
    y = ee.Image(y)
    num = x.gamma().multiply(y.gamma())
    den = x.add(y).gamma()
    return ee.Image(num.divide(den))


def incbetaNum(x, a, b):
    """
    Incomplete beta function (sans-gamma term), scalar version.
    @param {Number} x  upper limit of integration
    @param {Number} a  first shape parameter
    @param {Number} b  second shape parameter
    @returns {ee.Number} incomplete beta value
    """
    a = ee.Number(a)
    b = ee.Number(b)
    x = ee.Number(x)
    # Higher count yields more precise numbers but adds computation time
    count = 100
    v = ee.List.sequence(0, x, None, count)
    # Calculate midpoint of intervals
    s1 = ee.Array(v.slice(0, count - 1))
    s2 = ee.Array(v.slice(1, count))
    v = s1.add(s2).divide(2)
    # Get interval width
    dx = v.get([1]).subtract(v.get([0]))
    unit_array = ee.Array(ee.List.repeat(1, count - 1))
    _incbeta = v.pow(a.subtract(UNIT)).multiply(unit_array.subtract(v).pow(b.subtract(UNIT)))
    incbeta = _incbeta.reduce(ee.Reducer.sum(), [0]).multiply(dx).get([0])
    return incbeta


def incbeta(x, a, b):
    """
    Incomplete beta function, image version.
    @param {ee.Image} x  upper limit image
    @param {ee.Image} a  first shape parameter image
    @param {ee.Image} b  second shape parameter image
    @returns {ee.Image} incomplete beta image
    """
    a = ee.Image(a)
    b = ee.Image(b)
    x = ee.Image(x)
    n_int = ee.Number(100)
    dx = x.divide(n_int)
    seq = ee.Image.constant(ee.List.sequence(0, n_int))
    iv = dx.multiply(seq)
    # Calculate midpoint of intervals
    s1 = iv.select(ee.List.sequence(0, n_int.subtract(1)))
    s2 = iv.select(ee.List.sequence(1, n_int))
    iv = s1.add(s2).divide(2)
    unit_img = ee.Image(1)
    _incbeta = iv.pow(a.subtract(unit_img)).multiply(unit_img.subtract(iv).pow(b.subtract(unit_img)))
    incbeta_out = _incbeta.reduce(ee.Reducer.sum()).multiply(dx)
    return incbeta_out


def regincbetaNum(x, a, b):
    """
    Regularized incomplete beta function, scalar version.
    @param {Number} x  value to evaluate at
    @param {Number} a  first shape parameter
    @param {Number} b  second shape parameter
    @returns {ee.Number} regularized incomplete beta value
    """
    return ee.Number(incbetaNum(x, a, b)).divide(ee.Number(incbetaNum(1, a, b)))


def regincbeta(x, a, b):
    """
    Regularized incomplete beta function, image version.
    @param {ee.Image} x  value image
    @param {ee.Image} a  first shape parameter image
    @param {ee.Image} b  second shape parameter image
    @returns {ee.Image} regularized incomplete beta image
    """
    return ee.Image(incbeta(x, a, b)).divide(ee.Image(incbeta(1, a, b)))


def F_cdfNum(x, v1, v2):
    """
    CDF of F distribution, scalar version.
    @param {Number} x   value to evaluate at
    @param {Number} v1  numerator degrees of freedom
    @param {Number} v2  denominator degrees of freedom
    @returns {ee.Number} CDF value
    """
    x = ee.Number(x)
    v1 = ee.Number(v1)
    v2 = ee.Number(v2)
    k = v2.divide(v2.add(v1.multiply(x)))
    cdf = UNIT.subtract(regincbetaNum(k, v2.divide(2), v1.divide(2)))
    return cdf


def F_cdf(x, v1, v2):
    """
    CDF of F distribution, image version.
    @param {ee.Image} x   value image
    @param {ee.Image} v1  numerator degrees of freedom image
    @param {ee.Image} v2  denominator degrees of freedom image
    @returns {ee.Image} CDF image
    """
    x = ee.Image(x)
    v1 = ee.Image(v1)
    v2 = ee.Image(v2)
    k = v2.divide(v2.add(v1.multiply(x)))
    cdf = ee.Image(1).subtract(regincbeta(k, v2.divide(2), v1.divide(2)))
    return cdf


def F_pdf(x, df1, df2):
    """
    F distribution PDF using implementation from scipy.stats.f.
    @param {Number} x    value to evaluate at
    @param {Number} df1  numerator degrees of freedom
    @param {Number} df2  denominator degrees of freedom
    @returns {ee.Number} PDF value
    """
    x = ee.Number(x)
    df1 = ee.Number(df1)
    df2 = ee.Number(df2)
    num = (df2.pow(df2.divide(2))
              .multiply(df1.pow(df1.divide(2)))
              .multiply(x.pow(df1.divide(2).subtract(1))))
    den = (df2.add(df1.multiply(x))
              .pow(df1.add(df2).divide(2))
              .multiply(beta(df1.divide(2), df2.divide(2))))
    return ee.Number(num.divide(den))
