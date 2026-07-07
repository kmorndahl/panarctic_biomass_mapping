"""

DESCRIPTION: Functions to work with field biomass harvest data

FUNCTION LIST:
  buffer_func
  mean_sd_cv
  water_mask
  water_summary
  water_filt
  export_rep

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:

TO-DO:

"""

import ee

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# FUNCTION: buffer_func
# USE: Create a square bounding region around a plot/site point location
# PARAMETERS:
#  dataset_id = numerical dataset identifier
#  feat_col = feature collection of plots/sites corresponding to dataset_id
#  buffer_widths = dictionary of site specific buffer widths
# NOTES: Parameters are populated automatically within 1_proc_site_representativeness script
# AUTHOR: Melissa Rose, NAU
# LAST UPDATE: 7-19-2022
# TO-DO:

def buffer_func(dataset_id, feat_col):

    # Plot-level
    def buffer_plot(feature):
        return feature.buffer(30).bounds()  # Return bounds of buffered feature

    subset1 = feat_col.filter(ee.Filter.equals('coord_type', 'plot')).map(buffer_plot)

    # Site-level
    def buffer_site(feature):
        return feature.buffer(90).bounds()  # Return bounds of buffered feature

    subset2 = feat_col.filter(ee.Filter.equals('coord_type', 'site')).map(buffer_site)

    buffer_feat_col = subset1.merge(subset2)  # Combine plot and site level buffered points

    return buffer_feat_col


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# FUNCTION: mean_sd_cv
# USE: Calculate the mean, standard deviation and coefficient of variation near plots/sites for selected bands
# PARAMETERS:
#  dataset_id = numerical dataset identifier
#  feat_col = feature collection of buffered plot/site points corresponding to dataset_id
#  img = composite image, bands = list of bands over which to calculate CV
# NOTES: Parameters are populated automatically within 1_proc_site_representativeness script
# TO-DO: Enable band selection to change
# AUTHOR: Melissa Rose, NAU
# LAST UPDATE: 7-19-2022
# TO-DO:

def mean_sd_cv(dataset_id, feat_col, img, bands):

    # Specify and combine reducers
    reducers = ee.Reducer.mean().combine(
        reducer2=ee.Reducer.stdDev(),
        sharedInputs=True
    )

    # Get mean and std over buffered plot/site points
    dataset_summary = img.select(bands).reduceRegions(
        collection=feat_col,
        reducer=reducers,
        scale=10,
    )

    # Map over dataset summary to get coefficient of variation
    def calc_cv(feature):
        ndvi_mean = ee.Number(feature.get('NDVI_mean'))
        ndvi_sd = ee.Number(feature.get('NDVI_stdDev'))
        ndvi_coeff_var = ee.Algorithms.If(
            ee.Algorithms.IsEqual(ndvi_mean, None),   # Check if mean is null
            None,                                     # If null report null
            ee.Number(ndvi_sd.divide(ndvi_mean).multiply(100)).float()  # If not null, calculate CV
        )

        ndmi_mean = ee.Number(feature.get('NDMI_mean'))
        ndmi_sd = ee.Number(feature.get('NDMI_stdDev'))
        ndmi_coeff_var = ee.Algorithms.If(
            ee.Algorithms.IsEqual(ndmi_mean, None),   # Check if mean is null
            None,                                     # If null report null
            ee.Number(ndmi_sd.divide(ndmi_mean).multiply(100)).float()  # If not null, calculate CV
        )

        nbr_mean = ee.Number(feature.get('NBR_mean'))
        nbr_sd = ee.Number(feature.get('NBR_stdDev'))
        nbr_coeff_var = ee.Algorithms.If(
            ee.Algorithms.IsEqual(nbr_mean, None),    # Check if mean is null
            None,                                     # If null report null
            ee.Number(nbr_sd.divide(nbr_mean).multiply(100)).float()   # If not null, calculate CV
        )

        return (feature
                .set('ndvi_cv', ndvi_coeff_var)
                .set('ndmi_cv', ndmi_coeff_var)
                .set('nbr_cv', nbr_coeff_var))

    summary = dataset_summary.map(calc_cv)

    return summary


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# FUNCTION: water_mask
# USE: Create a binary layer classifying pixels with and without water, snow/ice and wetland
# PARAMETERS:
#  dataset_id = numerical dataset identifier
#  feat_col = feature collection of buffered plot/site points (with coefficient of variation data) corresponding to dataset_id
#  img_ref = reference image collection to define projection for land cover data
# NOTES: Parameters are populated automatically within 1_proc_site_representativeness script
# AUTHOR: Melissa Rose, NAU
# LAST UPDATE: 9-13-2022
# TO-DO:

def water_mask(dataset_id, feat_col, ic_ref):

    # Create a bounding box for entire study region
    bbox = ee.Geometry(feat_col.geometry().bounds())

    # Import Sentinel-derived landcover classification (10m-resolution)
    landcov = ee.ImageCollection('ESA/WorldCover/v100').first()
    landcov = ee.Image(landcov).clip(bbox)

    # Add classified image for study region to map

    # Create binary mask for 'Snow and ice' (70) class
    snow_ice = landcov.eq(70).rename('Snow and ice')

    # Create binary mask for 'Open water' (80) class
    water = landcov.eq(80).rename('Open water')

    # Create binary mask for 'Herbaceous wetland' (90) class
    wetland = landcov.eq(90).rename('Herbaceous wetland')

    # Combine all three masks into 3-band mask image
    new_bands = ee.Image([water, wetland])
    water_mask_img = snow_ice.addBands(new_bands)

    return water_mask_img


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# FUNCTION: water_summary
# USE: Summarize the presence of water pixels within plot/site buffered areas - # of pixels and fraction of pixels
# PARAMETERS:
#  dataset_id = numerical dataset identifier
#  feat_col = feature collection of buffered plot/site points (with coefficient of variation data) corresponding to dataset_id
#  img = multi-band binary water mask image
# NOTES: Parameters are populated automatically within 1_proc_site_representativeness script
# AUTHOR: Melissa Rose, NAU
# LAST UPDATE: 9-13-2022
# TO-DO:

def water_summary(dataset_id, feat_col, img):

    # Count total number of pixels in each plot/site
    pixel_count = img.select(0).reduceRegions(
        reducer=ee.Reducer.count(),
        collection=feat_col,
        scale=10,
    )

    # Create single mask for all bands
    # Output is a multi-band image with masks for water, snow/ice and wetland
    # Converts from binary 1/0 mask to 1s and NAs so that the count reducer will work
    multi_band_mask = img.eq(1)

    # Count number of pixels in each band for each plot/site
    # i.e. for each plot/site count the number of water, snow/ice and wetland pixels
    water_count = img.mask(multi_band_mask).reduceRegions(
        reducer=ee.Reducer.count(),
        collection=pixel_count,
        scale=10,
    )

    # Calculate the fraction of pixels in each band for each site/plot
    def calc_fraction(feature):
        snow_pixel = ee.Number(feature.get('Snow and ice'))
        water_pixel = ee.Number(feature.get('Open water'))
        wetland_pixel = ee.Number(feature.get('Herbaceous wetland'))
        total_pixel = ee.Number(feature.get('count'))

        snow_fraction = ee.Number(snow_pixel.divide(total_pixel).float())
        water_fraction_val = ee.Number(water_pixel.divide(total_pixel).float())
        wetland_fraction = ee.Number(wetland_pixel.divide(total_pixel).float())

        water_feat = feature.set(
            'total_pixel_count', total_pixel,
            'snow_pixel_count', feature.get('Snow and ice'),       # renaming each band's count output
            'water_pixel_count', feature.get('Open water'),
            'wetland_pixel_count', feature.get('Herbaceous wetland'),
            'snow_fraction', snow_fraction,                        # creating new property for each band's relative percentage
            'water_fraction', water_fraction_val,
            'wetland_fraction', wetland_fraction
        )
        water_feat_properties = water_feat.propertyNames()  # Get property names
        keep_properties = water_feat_properties.filter(
            ee.Filter.inList('item', ['Snow and ice', 'Open water', 'Herbaceous wetland', 'count']).Not()
        )  # Get list of properties to keep

        return water_feat.select(keep_properties)  # Remove unwanted properties

    water_fraction = water_count.map(calc_fraction)

    return water_fraction


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# FUNCTION: water_filt
# USE: Summarize the presence of water pixels within plot/site buffered areas - # of pixels and fraction of pixels
# PARAMETERS:
#  dataset_id = numerical dataset identifier
#  feat_col = feature collection of buffered plot/site points (with coefficient of variation data) corresponding to dataset_id
#  value = threshold of water pixels to use for filtering representative vs. non-representative plots/sites (unit is # of pixels)
#  folder = Google Drive folder to export to
# NOTES: Parameters are populated automatically within 1_proc_site_representativeness script
# AUTHOR: Melissa Rose, NAU
# LAST UPDATE: 9-13-2022
# TO-DO:

def water_filt(dataset_id, feat_col, value, folder):

    # Assign each plot/site as 'representative' or 'non-representative'
    def assign_rep(feature):
        water_count = ee.Number(feature.get('snow_pixel_count')).add(
            ee.Number(feature.get('water_pixel_count'))
        )  # Get total number of water and snow/ice pixels
        water_rep = ee.Algorithms.If(
            water_count.lte(ee.Number(value)),
            'representative',
            'non-representative'
        )  # If total number of water and snow/ice pixels is less than user specified threshold assign 'representative' otherwise assign 'non-representative'
        return feature.set('representativeness', water_rep)  # Create and populate 'representativeness' property

    site_rep = feat_col.map(assign_rep)

    # Visualize
    # subset_rep = site_rep.filter(ee.Filter.equals('representativeness', 'representative'))
    # subset_nonrep = site_rep.filter(ee.Filter.equals('representativeness', 'non-representative'))
    # empty = ee.Image().byte()
    # site_rep_outline = empty.paint(featureCollection=subset_rep, color=1, width=1)
    # site_nonrep_outline = empty.paint(featureCollection=subset_nonrep, color=1, width=1)

    return site_rep


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# FUNCTION: export_rep
# USE: Export the representativeness results
# PARAMETERS:
#  dataset_id = numerical dataset identifier
#  feat_col = feature collection of buffered plot/site points (with coefficient of variation data) corresponding to dataset_id
#  folder = Google Drive folder to export to
# NOTES: Parameters are populated automatically within 1_proc_site_representativeness script
# AUTHOR: Katie Orndahl, NAU
# LAST UPDATE: 3-29-2023
# TO-DO:

def export_rep(dataset_id, feat_col, folder):

    # Export as .csv
    ee.batch.Export.table.toDrive(
        collection=feat_col,
        description=str(dataset_id) + '_gee_representativeness',
        folder=folder,
        fileFormat='CSV'
    ).start()

    # Export as .shp
    ee.batch.Export.table.toDrive(
        collection=feat_col,
        description=str(dataset_id) + '_gee_representativeness',
        folder=folder,
        fileFormat='SHP'
    ).start()

    return feat_col
