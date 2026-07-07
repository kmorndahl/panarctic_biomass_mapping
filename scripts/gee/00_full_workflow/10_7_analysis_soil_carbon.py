"""

DESCRIPTION: Calculate soil carbon across the Pan Arctic

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

# Reprojection options
crs = params.crs
canopy_height_threshold = params.canopy_height_threshold
scale = params.scale

# 1.1 ----- READ IN DATA -----

total_biomass = ee.Image('projects/arctic-biomass-mapping/assets/modeled_final/total_biomass_2020_j_index_v20240514').select('total_biomass')
canopy_height = ee.ImageCollection("projects/meta-forest-monitoring-okw37/assets/CanopyHeight").mosaic()
zones_fc = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/zones').select('FIRST_zone', 'dsl')
soil_carbon = ee.Image("projects/soilgrids-isric/ocs_mean").multiply(100)  # Convert from t/ha to g/m2

# =========================
# 2. TIDY DATA ============
# =========================

# 2.0 ----- SOIL CARBON -----

# Area image
# Area will be calculated at whatever scale is specified for reduceRegion/export
area = ee.Image.pixelArea()

# Multiply by pixel area to get total grams in each pixel
soil_carbon_g = soil_carbon.multiply(area).regexpRename('$', '_g')

# Mask by canopy height
canopy_height = canopy_height.unmask(0)  # Fill currently masked areas with low value to ensure we include Greenland, barrier islands etc.
canopy_height_mask = canopy_height.lt(canopy_height_threshold).selfMask()
soil_carbon_g = soil_carbon_g.updateMask(canopy_height_mask)

# Mask permanent ice and water from soil carbon map to match masking in our maps
mask_soil_carbon = total_biomass.mask()
soil_carbon_g = soil_carbon_g.updateMask(mask_soil_carbon)

# 2.1 ----- BIOCLIMATE ZONES -----

# Tidy zones
zones_img = zones_fc.reduceToImage(['dsl'], ee.Reducer.min()).byte()

# Create arctic image
arctic_img = zones_img.selfMask()

# Get zone information
zone_codes = zones_fc.aggregate_array('dsl').add(4)
zone_dict = ee.Dictionary({'1': 'High Arctic', '2': 'Low Arctic', '3': 'Oro Arctic', '4': 'Pan Arctic'})

# =========================
# 3. SUMMARIZE ============
# =========================

def summarize_zone(zone_code):
    # Get zone name
    zone_name = zone_dict.get(ee.Number(zone_code).format('%.0f'))

    # Mask soil carbon to current code only
    soil_carbon_g_zone_code = ee.Image(ee.Algorithms.If(
        ee.Number(zone_code).eq(4),
        soil_carbon_g.updateMask(arctic_img),
        soil_carbon_g.updateMask(zones_img.eq(ee.Image.constant(zone_code)).selfMask())
    ))

    # Tidy mask
    soil_carbon_g_zone_code = soil_carbon_g_zone_code.mask(soil_carbon_g_zone_code.mask().gt(0))

    # Reduce over masked biomass image
    soil_carbon_g_zone_code_reduced = soil_carbon_g_zone_code.reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=ee.Geometry.Polygon(params.export_region_coords, None, False),
        scale=soil_carbon.projection().nominalScale(),
        crs=soil_carbon.projection(),
        maxPixels=1e13,
        tileScale=16  # Default: 1
    )

    # Convert to feature
    feat = ee.Feature(None, soil_carbon_g_zone_code_reduced).set('dsl', zone_code).set('FIRST_zone', zone_name)

    return feat

soil_carbon_zone_summary = list(map(summarize_zone, zone_codes.getInfo()))

# =========================
# 4. EXPORT ===============
# =========================

task = ee.batch.Export.table.toDrive(
    collection=ee.FeatureCollection(soil_carbon_zone_summary).select(['.*'], None, False),
    description='soil_carbon_summary_mask' + str(canopy_height_threshold) + '_' + version,
    folder='biomass_summaries'
)
task.start()

raise Exception('stop')
