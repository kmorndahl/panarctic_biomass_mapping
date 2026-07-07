"""

DESCRIPTION: Export final map mask to drive for use in calculating masked area

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

# 1.0 ----- PARAMETERS -----

# Model options
version = 'v20240514'

# Mask options
canopy_height_threshold = params.canopy_height_threshold

# Reprojection options
crs = params.crs
scale = params.scale  # Scale for reprojecting predictors

# 1.1 ----- READ IN DATA -----

total_biomass = ee.Image('projects/arctic-biomass-mapping/assets/modeled_final/total_biomass_2020_final_' + version).select('total_biomass')
canopy_height = ee.ImageCollection("projects/meta-forest-monitoring-okw37/assets/CanopyHeight").mosaic()

# ===============================
# 2. APPLY CANOPY HEIGHT MASK ===
# ===============================

canopy_height = canopy_height.unmask(0)  # Fill currently masked areas with low value to ensure we include Greeland, barrier islands etc.
canopy_height_mask = canopy_height.lt(canopy_height_threshold).selfMask()
total_biomass_mask = total_biomass.updateMask(canopy_height_mask).multiply(0).add(1)


# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.image.toDrive(
    image=total_biomass_mask.byte(),  # Convert to byte for easier storage
    description='final_biomass_mask_trees',
    folder='biomass_summaries',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    maxPixels=1e13
)
task.start()

raise Exception('stop')