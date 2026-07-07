"""

DESCRIPTION: Create permafrost predictor

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:
- Tips for exporting global asset:
  - https://groups.google.com/g/google-earth-engine-developers/c/55btJBTGQoI/m/XUZWpLG1EAAJ
  - https://groups.google.com/g/google-earth-engine-developers/c/SvCgqmXbMS0/m/KcoziJugCwAJ

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

# 1.0 ----- READ IN DATA -----

# Permafrost Zonation Index (Gruber 2012)
permafrostMetrics = ee.Image('projects/foreststructure/misc/global_permafrost_zonation_index_gruber2012')

# 1.1 ----- PARAMETERS -----

# Projection (derived from input data)
proj = permafrostMetrics.projection().getInfo()
crs = proj['crs']
transform = proj['transform']

# =========================
# 2. ANALYSIS =============
# =========================

# Original was 0.008333 degrees, caused rounding error after 43200 pixels
permafrostMetrics_shifted = permafrostMetrics.changeProj(
    permafrostMetrics.projection(),
    ee.Projection('EPSG:4326', [0.0083333333333333, 0, -180, 0, -0.0083333333333333, 89.994])
)

# Calculate 20 km focal mean, leaves sizeable gap
permafrostMetrics_shifted_fill = permafrostMetrics_shifted.focal_mean(20000, 'circle', 'meters')

# Create geometry capturing gap
geom_gap = ee.FeatureCollection(ee.Geometry.BBox(-179, -89.99, -177, 89.99))

# Clip original image to extent of gap
permafrostMetrics_shifted_gap = permafrostMetrics_shifted.clip(geom_gap)

# Calculate 5 km focal mean on gap image
permafrostMetrics_shifted_gap_filled = permafrostMetrics_shifted_gap.focal_mean(5000, 'circle', 'meters')
permafrostMetrics_shifted_gap_filled = permafrostMetrics_shifted_gap_filled.updateMask(
    permafrostMetrics_shifted_gap_filled.gt(0.5)
)

# Combine focal mean images
permafrostMetrics_shifted_fill_all = ee.ImageCollection([
    permafrostMetrics_shifted_fill,
    permafrostMetrics_shifted_gap_filled
]).mosaic()

# Add back in actual data
permafrostMetrics_final = ee.ImageCollection([
    permafrostMetrics_shifted_fill_all,
    permafrostMetrics_shifted
]).mosaic().unmask(0)

# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.image.toAsset(
    image=permafrostMetrics_final,
    description='permafrostMetrics_final',
    assetId='projects/arctic-biomass-mapping/assets/land_cover/global_permafrost_zonation_index_gruber2012_filled_20km_v3',
    dimensions='43200x18000',
    crs=crs,
    crsTransform=[0.0083333333333333, 0, -180, 0, -0.0083333333333333, 89.994],
    maxPixels=1e12
)
task.start()

raise Exception('stop')
