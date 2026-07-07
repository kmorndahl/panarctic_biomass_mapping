"""

DESCRIPTION: Create world terrestrial ecosystem predictor

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

# 1.0 ----- READ IN DATA ----

# Get World Terrestrial Ecosystems data
# https://www.sciencebase.gov/catalog/item/6296791ed34ec53d276bb293
# https://rmgsc.cr.usgs.gov/outgoing/ecosystems/Global/
wte = ee.Image('projects/arctic-biomass-mapping/assets/predictors/wte_2020').rename('world_terrestrial_ecosystems')

# Get Greenland data
greenland_water = ee.Image("OSU/GIMP/2000_ICE_OCEAN_MASK").select('ocean_mask')

# 1.1 ----- PARAMETERS ----

# Projection (derived from input data)
proj = wte.projection()
crs = proj.crs()
scale = proj.nominalScale()
transform = proj.getInfo()['transform']
transform_scale = transform[0]

# =========================
# 2. TIDY =================
# =========================

# 2.0 ----- FILL GAPS - FOCAL MODE ----

# Fill some gaps to avoid NAs in cal/val data
wte_fill_gaps = wte.focalMode(
    radius=scale.multiply(3).divide(2),
    kernelType='square',
    units='meters',
    iterations=3
)
wte = wte_fill_gaps.blend(wte)
print('World Terrestrial Ecosystems - Filled', wte.getInfo())

# 2.1 ----- FILL GREENLAND GAPS ----

# Tidy greenland data
greenland_land = greenland_water.eq(0).selfMask()

# Split Greenland
# S will be assigned WTE class 2 = "Polar Moist Snow and Ice on Plains" to fill large ice gap
# N will be assigned WTE class 1 = "Polar Moist Sparsely or Non-vegetated on Plains" to account for rocky coastlines
n_greenland_poly = ee.Geometry.Polygon([-89, 88, 0, 88, 7.6, 88, 7.6, 78.7, 0, 78.7, -89, 78.7], None, False)
s_greenland_poly = ee.Geometry.Polygon([-89, 78.7, 0, 78.7, 7.6, 78.7, 7.6, 58.8, 0, 58.8, -89, 58.8], None, False)
n_greenland = greenland_land.clip(n_greenland_poly)
s_greenland = greenland_land.clip(s_greenland_poly).multiply(2)
greenland_land = n_greenland.blend(s_greenland)

# Reproject greenland to match WTE
greenland_land = greenland_land.reproject(
    crs='EPSG:4326',
    crsTransform=[transform_scale, 0, -180, 0, -transform_scale, 90]
)

# Fill Greenland gaps
wte = greenland_land.blend(wte).rename('world_terrestrial_ecosystems')
wte = wte.updateMask(wte.mask().gt(0))  # Remove mask percentage
print('World Terrestrial Ecosystems - Greenland Filled', wte.getInfo())

# 2.2 ----- EXPORT ----

task = ee.batch.Export.image.toAsset(
    image=wte,
    description='wte_2020_filled',
    assetId='projects/arctic-biomass-mapping/assets/predictors/wte_2020_filled',
    pyramidingPolicy='mode',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    crs='EPSG:4326',
    crsTransform=[transform_scale, 0, -180, 0, -transform_scale, 90],
    maxPixels=1e12
)
task.start()

# =========================
# 3. ENCODE ===============
# =========================

# 3.0 ----- ENCODE ----

# Remove mask opacity
wte = wte.updateMask(wte.mask().gt(0))

# One-hot encode
hist = wte.reduceRegion(
    reducer=ee.Reducer.frequencyHistogram(),
    geometry=ee.Geometry.Polygon(params.export_region_coords, None, False),
    maxPixels=1e13
)
# Remove all the unnecessary reducer output structure and make a list of values
values = ee.List(ee.Dictionary(hist.get(wte.bandNames().get(0))).keys().map(ee.Number.parse)).getInfo()

def make_wte_band(value):
    return wte.eq(ee.Image.constant(value)).rename('world_terrestrial_ecosystems_X' + str(int(value)))

wte = (ee.ImageCollection(list(map(make_wte_band, values)))
       .toBands()
       .regexpRename('^[^_]*_', '')  # Remove prefixes added by 'toBands()'
       .uint8())
print('World Terrestrial Ecosystems - Encoded', wte.getInfo())

# 3.1 ----- EXPORT ----

# Export
task = ee.batch.Export.image.toAsset(
    image=wte,
    description='wte_2020_filled_encoded',
    assetId='projects/arctic-biomass-mapping/assets/predictors/wte_2020_filled_encoded',
    pyramidingPolicy='mode',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    crs='EPSG:4326',
    crsTransform=[transform_scale, 0, -180, 0, -transform_scale, 90],
    maxPixels=1e12
)
task.start()

raise Exception('stop')


