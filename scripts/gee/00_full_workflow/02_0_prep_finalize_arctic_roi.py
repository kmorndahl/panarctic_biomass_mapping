"""

DESCRIPTION: Finalize Arctic ROI by manually removing lingering holes and crumbs

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:
- Most processing of the Arctic ROI is done in R
- This script performs the final, manual tidying steps

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

arctic_oroarctic_coast_buffer_laea = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/arctic_oroarctic_coast_buffer_laea')  # Buffered collection with some unwanted holes and crumbs
arctic_oroarctic_fill_holes = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/arctic_oroarctic_fill_holes')  # User created feature collection to fill holes
arctic_oroarctic_remove_crumbs = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/arctic_oroarctic_remove_crumbs')  # User created feature collection to remove crumbs

# 1.1 ----- PARAMETERS -----

# Reprojection options
crs = params.crs

# =========================
# 2. TIDY =================
# =========================

# 2.0 ----- FILL IN HOLES -----

arctic_tidy = arctic_oroarctic_coast_buffer_laea.merge(ee.FeatureCollection(arctic_oroarctic_fill_holes)).union(30)

# 2.1 ----- REMOVE CRUMBS -----

# Separate crumbs
lac_brochet = arctic_oroarctic_remove_crumbs.filter(ee.Filter.eq('system:index', '00000000000000000000')).first()
swe = arctic_oroarctic_remove_crumbs.filter(ee.Filter.eq('system:index', '00000000000000000001')).first()
flin_flon = arctic_oroarctic_remove_crumbs.filter(ee.Filter.eq('system:index', '00000000000000000002')).first()

# Remove crumbs
arctic_rm_lac = ee.Feature(arctic_tidy.first()).difference(ee.Feature(lac_brochet), maxError=30, proj=crs)
arctic_rm_swe = ee.Feature(arctic_rm_lac).difference(ee.Feature(swe), maxError=30, proj=crs)
arctic_final = ee.Feature(arctic_rm_swe).difference(ee.Feature(flin_flon), maxError=30, proj=crs)

# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.table.toAsset(
    collection=ee.FeatureCollection([arctic_final]),
    description='arctic_oroarctic_coast_buffer_laea_final',
    assetId='projects/arctic-biomass-mapping/assets/ROIs/arctic_oroarctic_coast_buffer_laea_final'
)
task.start()

raise Exception('stop')
