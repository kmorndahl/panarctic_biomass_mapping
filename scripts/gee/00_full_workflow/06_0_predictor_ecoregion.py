"""

DESCRIPTION: Create ecoregion predictor

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:

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

# 1.0 ----- PARAMETERS ----

# Reprojection options
crs = params.crs
scale = params.scale  # Scale for reprojecting predictors

# 1.1 ----- READ IN DATA ----

ecoregion = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/ecoregions_resolve_2017_45n')

# 1.2 ----- TIDY DATA ----

# Tidy ecoregion names
def tidy_ecoregion_name(feat):
    return feat.set({'ECO_NAME': ee.String(feat.get('ECO_NAME'))
                                   .replace(' ', '', 'g')
                                   .replace('[.]', '', 'g')
                                   .replace('-', '', 'g')
    })

ecoregion = ecoregion.map(tidy_ecoregion_name)

print('Ecoregion feature collection:', ecoregion.getInfo())

# =========================
# 2. ANALYSIS =============
# =========================

# 2.0 ----- GET ECOREGION NAMES AND IDs ----

# Get ecoregion names and IDs from featureCollection
ecoregion_names = ecoregion.aggregate_array('ECO_NAME')
ecoregion_nums = ecoregion.aggregate_array('ECO_ID')
num_ecoregion = 846  # https://www.arcgis.com/home/item.html?id=37ea320eebb647c6838c23f72abae5ef

# 2.1 ----- CONVERT ECOREGION VECTOR FILE TO IMAGE ----

ecoregion_img = ecoregion.reduceToImage(
    properties=['ECO_ID'],
    reducer=ee.Reducer.first()
)

# All masked areas are converted to 847 to represent 'ecoregion_other'
ecoregion_img = ecoregion_img.unmask(ee.Number(num_ecoregion).add(1))

# Add 'other' class
ecoregion_nums = ecoregion_nums.cat([ee.Number(num_ecoregion).add(1)])  # Create list of ecoregion numbers, including unique number for 'ecoregion_other'
ecoregion_names = ecoregion_names.cat(['other']).map(lambda name: ee.String('ecoregion_').cat(name))  # Create list of ecoregion names, including unique name for 'ecoregion_other'
print('Ecoregion names:', ecoregion_names.getInfo())
print('Ecoregion numbers:', ecoregion_nums.getInfo())

# 2.2 ----- CONVERT SINGLE BAND ECOREGION IMAGE INTO MULTI-BAND BINARY IMAGE ----

def make_ecoregion_band(num):
    return ecoregion_img.eq(ee.Image.constant(num)).unmask(0)

ecoregion_final = ee.ImageCollection(ecoregion_nums.map(make_ecoregion_band))

ecoregion_final = (ecoregion_final
    .toBands()
    .rename(ecoregion_names)
    .setDefaultProjection(crs=crs, scale=scale)
    .uint8())

print(ecoregion_final.getInfo())

# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.image.toAsset(
    image=ecoregion_final,
    description='ecoregions_img',
    assetId='projects/arctic-biomass-mapping/assets/predictors/ecoregions_img',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    maxPixels=1e13,
    pyramidingPolicy={'.default': 'mode'}
)
task.start()

raise Exception('stop')
