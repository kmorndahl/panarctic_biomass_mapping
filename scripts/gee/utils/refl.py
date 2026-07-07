import ee

def fitMapT(doyRaster, year, segment, tag, extrapolateMaxDays, shown, Map=None, visParams=None):
    """
    Model reflectance over a DOY image, and map.
    """
    # Combine year and DOY to specify exact date (fractional year)
    t = ee.Image(year).add(doyRaster.divide(365.25))
    
    # Python API uses keyword arguments for dictionary-style parameters
    fit = segment.slice(t=t, harmonics=3, extrapolateMaxDays=extrapolateMaxDays)
    
    if Map:
        Map.addLayer(fit, visParams, f'fit {tag}', shown)
        
    return fit

def fitT(doyRaster, year, segment, extrapolateMaxDays):
    """
    Model reflectance over a DOY image.
    """
    t = ee.Image(year).add(doyRaster.divide(365.25))
    fit = segment.slice(t=t, harmonics=3, extrapolateMaxDays=extrapolateMaxDays)
    return fit

def NDVIslope(segment, extrapolateMaxDays):
    """
    Calculate NDVI slope over the current CCDC segment.
    """
    segment_image = segment.segment_image
    segment_start = segment_image.select('tStart')
    segment_end = segment_image.select('tEnd')
    
    # Calculate length of segment in days
    segment_length = segment_end.subtract(segment_start).multiply(365.25).rename('tLength_days')
    
    # Model reflectance for segment start and end
    segment_start_fit = segment.slice(t=segment_start, harmonics=3, extrapolateMaxDays=extrapolateMaxDays)
    segment_end_fit = segment.slice(t=segment_end, harmonics=3, extrapolateMaxDays=extrapolateMaxDays)
    
    # Calculate change in NDVI per day, scaled by 10,000
    slope = (segment_end_fit.select('NDVI')
             .subtract(segment_start_fit.select('NDVI'))
             .divide(segment_length)
             .multiply(10000)
             .rename('ccdc_NDVIslope10000'))
    
    return slope

def calculateTexture(img, texture_radius, texture_bands, texture_metrics, crs, scale):
    """
    Calculate NDVI texture over the current CCDC segment.
    """
    kernel = ee.Kernel.square(radius=texture_radius)
    
    txt = (img.select(texture_bands)
           .setDefaultProjection(crs=crs, scale=scale)
           .unitScale(-10000, 10000)
           .multiply(255)
           .toByte()
           .glcmTexture(kernel=kernel, average=True)
           .select(texture_metrics)
           .regexpRename('spectral_', '')
           .regexpRename(ee.String(texture_metrics[0]).replace('[.][*]', ''), '')
           .regexpRename('^', 'texture_'))
    
    return txt

def calculateTimeSinceBreak(segment, date):
    """
    Calculate time since segment start.
    Note: Requires a date conversion logic. 
    """
    # If date is an ee.Date, this converts it to fractional year
    # Replacing the 'fun_misc.dateConversion.toT' logic:
    year = date.get('year')
    doy = date.getRelative('day', 'year')
    dateT = year.add(doy.divide(365.25))
    
    segment_start = segment.toImage().select('tStart')
    time_since_break = ee.Image.constant(dateT).subtract(segment_start)
    
    return time_since_break