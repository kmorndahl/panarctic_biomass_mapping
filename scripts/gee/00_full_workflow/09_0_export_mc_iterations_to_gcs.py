"""

DESCRIPTION: Export Monte Carlo iterations to Google Cloud Storage

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:

IMPORTANT:
- For best results, GCS buckets should be in US-CENTRAL1 (GEE servers are in Iowa)
- https://groups.google.com/g/google-earth-engine-developers/c/3MiE78ad8Aw/m/6F0Wi0XaAAAJ

TO-DO:

"""

import ee
from utils import params

try:
    ee.Initialize(project=params.ee_project)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=params.ee_project)

# =========================
# 1. SET-UP ===============
# =========================

# 1.0 ----- PARAMETERS -----

# Model options
response_type = 'continuous'
ds_type = 'woody'
year = 2020
version = 'v20240514'

# Reprojection options
crs = params.crs
scale = params.scale

# Export format options
no_data_val = -1 if response_type == 'binary' else 65535 # Choose -1 for binary, 65535 for continuous
file_dim = 131072 if response_type == 'binary' else 89600 # Choose 131072 for binary, 89600 for continuous

# 1.1 ----- READ IN DATA -----

ic = ee.ImageCollection('projects/arctic-biomass-mapping/assets/modeled_' + response_type + '/model_' + response_type + '_' + str(year) + '_' + ds_type + '_' + version + '_mc')

# =========================
# 2. EXPORT ===============
# =========================

# 2.0 ----- SET UP LIST OF MONTE CARLO ITERATIONS (AS STRINGS) -----

mc_iterations = ee.List.sequence(1, 100, 1).map(lambda item: ee.Number(item).format('%.0f')).getInfo()

# 2.1 ----- LOOP IMAGES IN IMAGE COLLECTION -----

for mc in mc_iterations:

    # Get image
    img = ic.filter(ee.Filter.eq('system:index', mc)).first()

    # Set no data values
    img = img.unmask(no_data_val)  # Set masked areas so they aren't automatically set to 0

    # Export
    task = ee.batch.Export.image.toCloudStorage(
        image=img,
        description=mc,
        bucket='mc_iterations',
        fileNamePrefix=str(year) + '/' + response_type + '/' + ds_type + '/mc' + mc + '_',
        region=ee.Geometry.Polygon(params.export_region_coords, None, False),
        scale=scale,
        crs=crs,
        maxPixels=1e13,
        fileDimensions=file_dim,  # https://gis.stackexchange.com/questions/359974/calculating-max-shardsize-and-filedimensions-to-avoid-tiling-large-raster-export/360062#360062
        formatOptions={  # https://developers.google.com/earth-engine/guides/exporting_images
            'cloudOptimized': True,
            'noData': no_data_val
        }
    )
    task.start()

raise Exception('stop')