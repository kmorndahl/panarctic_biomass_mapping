"""

DESCRIPTION: Average biomass across growing degree day bins

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

# Climate analysis options
gdd_base = 0
gdd_bin_width = 10 if gdd_base == 0 else 5  # For base 0, choose 10; for base 5, choose 2

# Mask options
canopy_height_threshold = params.canopy_height_threshold

# Output options
output_folder = 'biomass_mean'

# 1.1 ----- READ IN DATA -----

# CHELSA growing degree days:
#  - units: degrees celsius
#  - description: heat sum of all days above the 0 degrees C temperature accumulated over 1 year.
#  - scale: 0.1
gdd = ee.Image("projects/arctic-biomass-mapping/assets/climate_analysis/CHELSA_gdd" + str(gdd_base) + "_1981_2010_V_2_1").multiply(0.1).rename('gdd').round().unmask(0, False)

total_biomass = ee.Image("projects/arctic-biomass-mapping/assets/modeled_final/total_biomass_2020_j_index_v20240514")
woody_biomass = ee.Image("projects/arctic-biomass-mapping/assets/modeled_final/woody_biomass_2020_j_index_v20240514")
canopy_height = ee.ImageCollection("projects/meta-forest-monitoring-okw37/assets/CanopyHeight").mosaic()

# =========================
# 2. ANALYSIS =============
# =========================

# 2.0 ----- CALCULATE WOODY PLANT DOMINANCE -----

woody_percent = woody_biomass.divide(total_biomass).multiply(100).select('woody_biomass').regexpRename('biomass', 'percent')  # Where plant and woody plant biomass are both zero, woody percent is zero

# 2.1 ----- SET UP GDD DATA -----

# Get max GDD -----

gdd_max = gdd.updateMask(total_biomass.select('total_biomass')).reduceRegion(
    reducer=ee.Reducer.max(),
    geometry=ee.Geometry.Polygon(params.export_region_coords, None, False),
    scale=gdd.projection().nominalScale(),
    crs=gdd.projection().crs(),
    maxPixels=1e13
)

gdd_max = ee.Number(ee.Dictionary(gdd_max).get('gdd'))
gdd_max = gdd_max.divide(gdd_bin_width).round().multiply(gdd_bin_width).getInfo()

# Bin GDD data -----

thresholds = ee.List.sequence(gdd_bin_width, gdd_max, gdd_bin_width).getInfo()  # Upper limit is maximum across Arctic domain
threshold_img = ee.Image(thresholds)
gdd_bins = gdd.gt(threshold_img).reduce('sum')  # Pixel values are bin numbers
gdd_bins = gdd_bins.add(1).multiply(gdd_bin_width).rename('gdd_bin')  # Convert pixel values to bin end value

# Mask GDD data -----

# Mask to Arctic domain
gdd_bins = gdd_bins.updateMask(total_biomass.select('total_biomass').mask())

# Mask by canopy height
canopy_height = canopy_height.unmask(0)  # Fill currently masked areas with low value to ensure we include Greeland, barrier islands etc.
canopy_height_mask = canopy_height.lt(canopy_height_threshold).selfMask()
gdd_bins = gdd_bins.updateMask(canopy_height_mask)

# 2.2 ----- SET UP BIOMASS DATA -----

# Mask by canopy height -----
canopy_height = canopy_height.unmask(0)  # Fill currently masked areas with low value to ensure we include Greeland, barrier islands etc.
canopy_height_mask = canopy_height.lt(canopy_height_threshold).selfMask()
total_biomass = total_biomass.updateMask(canopy_height_mask)
woody_biomass = woody_biomass.updateMask(canopy_height_mask)

# 2.3 ----- COMBINE BIOMASS DATA -----

biomass = total_biomass.addBands(woody_biomass).addBands(woody_percent)

# 2.4 ----- AVERAGE BIOMASS ACROSS BINS AND EXPORT -----

for bin_val in thresholds:

    # Restrict GDD bin image to current bin only
    gdd_bin = gdd_bins.eq(ee.Image.constant(bin_val)).selfMask()

    # Mask biomass to current GDD bin only
    biomass_gdd_bin = biomass.updateMask(gdd_bin)

    # Tidy mask
    biomass_gdd_bin = biomass_gdd_bin.mask(biomass_gdd_bin.mask().gt(0))

    # Total
    biomass_gdd_bin_reduced = biomass_gdd_bin.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=ee.Geometry.Polygon(params.export_region_coords, None, False),
        crs=gdd.projection(),
        tileScale=1,
        maxPixels=1e13
    )

    biomass_per_gdd_bin_fc = ee.FeatureCollection(ee.Feature(None, biomass_gdd_bin_reduced).set('gdd_bin', bin_val))

    task = ee.batch.Export.table.toDrive(
        collection=biomass_per_gdd_bin_fc,
        description='biomass_gdd' + str(gdd_base) + '_bin' + str(bin_val) + '_mean_1km',
        folder=output_folder
    )
    task.start()

raise Exception('stop')

