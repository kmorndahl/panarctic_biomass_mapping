"""

DESCRIPTION: Create masks to remove ocean and ice

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

# 1.0 ----- SET PARAMETERS -----

# JRC water mask options
occurrence_percent = 90  # Percent threshold to use for JRC water mask
buffer_dist = 0  # Buffer distance for JRC water mask

# Reprojection options
crs = params.crs
scale = params.scale

# Output options
out_version = 'v20240206'

# 1.1 ----- READ IN AND TIDY DATA -----

# JRC observation data
# Must use earlier version for differentiating ocean/land
#  - JRC/GSW1_0/Metadata outdated and inaccurate for Greenland - we will delineate Greenland ourselves
#  - JRC/GSW1_4/Metadata has observations even in the ocean
jrc_obs = ee.Image('JRC/GSW1_0/Metadata').select('total_obs')

# WorldCover data
lc = ee.ImageCollection('ESA/WorldCover/v200').first()

# Greenland ocean data
greenland_ocean = (ee.Image('OSU/GIMP/2000_ICE_OCEAN_MASK')
                   .select('ocean_mask')
                   .rename('water_mask')
                   .eq(0)
                   .unmask(value=0, sameFootprint=False))

# Greenland ice data
greenland_ice = (ee.Image('OSU/GIMP/2000_ICE_OCEAN_MASK')
                 .select('ice_mask')
                 .rename('water_mask')
                 .eq(1)
                 .unmask(value=0, sameFootprint=False))

# =========================
# 2. CREATE MASKS =========
# =========================

# 2.0 ----- DEEP OCEAN AND ICESHEET MASK -----

# Get JRC observation mask
jrc_mask = jrc_obs.mask().rename(['water_mask'])

# Get WorldCover ocean mask
lc_mask = (lc.mask()  # WorldCover data includes some nearshore ocean labeled water, rest of ocean is masked
              .eq(1)
              .rename(['water_mask']))

# Create deep ocean mask
ocean_mask = (jrc_mask.add(lc_mask)  # Removes some lingering unmasked discs
                      .add(lc_mask.clip(ee.Geometry.BBox(-180, 78, 180, 90)).unmask(value=0, sameFootprint=False))  # Necessary because JRC data only goes to 78 degrees north
                      .gte(2)
                      .add(greenland_ocean)  # Fix Greenland
                      .gte(1))

# Add Greenland icesheet mask
ocean_ice_mask = (ocean_mask.subtract(greenland_ice)
                             .gte(1))

# Export
task = ee.batch.Export.image.toAsset(
    image=ocean_ice_mask.toInt8(),
    description='ocean_ice_mask',
    assetId='projects/arctic-biomass-mapping/assets/ROIs/ocean_ice_mask_' + out_version,
    region=ee.Geometry.BBox(-180, 45, 180, 90),
    scale=scale,
    crs=crs,
    maxPixels=1e13
)
task.start()

raise Exception('stop')
