"""

DESCRIPTION: Create tree cover predictors

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

# 1.0 ----- READ IN DATA -----

tree_cover = ee.Image("UMD/hansen/global_forest_change_2022_v1_10")

# 1.1 ----- PARAMETERS -----

# Processing year
year = 2020  # Select year for predictions

# Projection (derived from input data)
proj = tree_cover.projection()
crs = proj.crs()
transform = proj.getInfo()['transform']

# 1.2 ----- FUNCTIONS -----

# FUNCTION: calculateTrees2020
# USE: Calculate tree percent cover in the year 2020
# AUTHOR: Kathleen Orndahl
# LAST UPDATE: 11-10-2024
def calculateTrees2020(img):

    lost_before_2020 = (img.select('trees_lossyear')
        .unmask(100)  # Areas where there was no tree loss are set to very high value so they aren't assigned a loss
        .gte(20))

    tree_cover_updated = img.select('trees_cover').multiply(lost_before_2020).unmask(0)

    tree_presence_updated = img.select('trees_presence').add(img.select('trees_gain')).gt(0)

    return tree_cover_updated.addBands(tree_presence_updated)

# =========================
# 2. ANALYSIS =============
# =========================

# Tidy
tree_cover = (tree_cover
    .select(['treecover2000', 'loss', 'gain', 'lossyear'], ['cover', 'loss', 'gain', 'lossyear'])  # Select relevant bands
    .regexpRename('^', 'trees_'))  # Add prefix

# Calculate tree presence/absence
tree_cover = tree_cover.addBands(
    tree_cover.select('trees_cover').gt(0).rename('trees_presence')
)  # Calculate presence/absence

# Calculate tree cover for specified year
tree_cover = ee.Image(
    ee.Algorithms.If(
        ee.Number(year).eq(2020),
        calculateTrees2020(tree_cover),
        tree_cover.select(['trees_cover', 'trees_presence'])
    )
)  # Get tree predictors based on year

# Unmask so that high Arctic and ocean are filled with zeros
tree_cover = tree_cover.unmask(0)

# Check
print('Tree cover:', tree_cover.getInfo())

# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.image.toAsset(
    image=tree_cover,
    description='tree_cover_' + str(year),
    assetId='projects/arctic-biomass-mapping/assets/predictors/tree_cover_' + str(year),
    pyramidingPolicy='mean',
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    crs='EPSG:4326',
    crsTransform=transform,
    maxPixels=1e12
)
task.start()

raise Exception('stop')