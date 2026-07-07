"""

DESCRIPTION: Average biomass and GDD across CAVM vegetation community types

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
gdd_bin_width = 10

# Mask options
canopy_height_threshold = params.canopy_height_threshold

# Output options
output_folder = 'gdd' + str(gdd_base)

# 1.1 ----- READ IN DATA -----

# CHELSA growing degree days:
#  - units: degrees celsius
#  - description: heat sum of all days above the 0 degrees C temperature accumulated over 1 year.
#  - scale: 0.1
gdd = ee.Image("projects/arctic-biomass-mapping/assets/climate_analysis/CHELSA_gdd" + str(gdd_base) + "_1981_2010_V_2_1").multiply(0.1).rename('gdd').round().unmask(0, False)

total_biomass = ee.Image("projects/arctic-biomass-mapping/assets/modeled_final/total_biomass_2020_j_index_v20240514")
woody_biomass = ee.Image("projects/arctic-biomass-mapping/assets/modeled_final/woody_biomass_2020_j_index_v20240514")
cavm = ee.Image('projects/arctic-biomass-mapping/assets/ROIs/cavm_raster')
cavm_legend = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/cavm_raster_legend')
canopy_height = ee.ImageCollection("projects/meta-forest-monitoring-okw37/assets/CanopyHeight").mosaic()

# =========================
# 2. ANALYSIS =============
# =========================

# 2.0 ----- CALCULATE WOODY PLANT DOMINANCE -----

woody_percent = woody_biomass.divide(total_biomass).multiply(100).select('woody_biomass').regexpRename('biomass', 'percent')  # New option 'set0' - don't mask, where plant and woody plant biomass are both zero, woody percent is zero

# 2.1 ----- SET UP GDD DATA -----

# Get max GDD -----

gdd_max = gdd.updateMask(woody_biomass.select('woody_biomass')).reduceRegion(
    reducer=ee.Reducer.max(),
    geometry=ee.Geometry.Polygon(params.export_region_coords, None, False),
    scale=gdd.projection().nominalScale(),
    crs=gdd.projection().crs(),
    maxPixels=1e13
)

gdd_max = ee.Number(ee.Dictionary(gdd_max).get('gdd'))
gdd_max = gdd_max.divide(gdd_bin_width).round().multiply(gdd_bin_width).getInfo()

# 2.2 ----- SET UP BIOMASS DATA -----

# Mask by canopy height -----
canopy_height = canopy_height.unmask(0)  # Fill currently masked areas with low value to ensure we include Greeland, barrier islands etc.
canopy_height_mask = canopy_height.lt(canopy_height_threshold).selfMask()
total_biomass = total_biomass.updateMask(canopy_height_mask)
woody_biomass = woody_biomass.updateMask(canopy_height_mask)

# 2.3 ----- COMBINE DATA -----

biomass_gdd = total_biomass.addBands(woody_biomass).addBands(woody_percent).addBands(gdd)

# 2.4 ----- AVERAGE BIOMASS AND GDD ACROSS CAVM VEGETATION COMMUNITY TYPES AND EXPORT -----

# Get CAVM codes
cavm_codes_fine = cavm_legend.aggregate_array('raster_code')

# Pre-fetch CAVM legend as a Python dict before the loop.
# Each feature in cavm_legend has 'raster_code' and 'short_description' properties; fetching
# all at once avoids a separate GEE filter expression per code inside the loop.
cavm_dict = {
    feat['properties']['raster_code']: feat['properties']['short_description']
    for feat in cavm_legend.getInfo()['features']
}

# Loop CAVM codes
biomass_gdd_cavm_fine = []
for code in cavm_codes_fine.getInfo():

    # Get vegetation description
    veg_desc = cavm_dict[code]

    # Restrict CAVM image to current code only
    cavm_code = cavm.eq(ee.Image.constant(code)).selfMask()

    # Mask biomass to current code only
    biomass_gdd_cavm_code = biomass_gdd.updateMask(cavm_code)

    # Tidy mask
    biomass_gdd_cavm_code = biomass_gdd_cavm_code.mask(biomass_gdd_cavm_code.mask().gt(0))

    # Reduce over masked biomass image
    # Reduces using CRS and resolution from gdd due to crs: gdd.projection()
    biomass_gdd_cavm_code_reduced = biomass_gdd_cavm_code.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=ee.Geometry.Polygon(params.export_region_coords, None, False),
        crs=gdd.projection(),
        maxPixels=1e13,
        tileScale=1)  # Default: 0.5

    # Convert to feature
    feat = ee.Feature(None, biomass_gdd_cavm_code_reduced).set('veg_code', code).set('veg_desc', veg_desc)

    biomass_gdd_cavm_fine.append(feat)

# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.table.toDrive(
    collection=ee.FeatureCollection(biomass_gdd_cavm_fine).select(['.*'], None, False),
    description='arctic_biomass_gdd' + str(gdd_base) + '_mean_1km_cavm_fine',
    folder=output_folder
)
task.start()

raise Exception('stop')
