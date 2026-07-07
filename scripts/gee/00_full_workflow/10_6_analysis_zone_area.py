"""

DESCRIPTION: Calculate area across bioclimate zones

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:
- Country FeatureCollections are too complex -- area must be calculated in R

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

# 1.0 ----- READ IN AND TIDY DATA -----

total_biomass = ee.Image('projects/arctic-biomass-mapping/assets/modeled_final/total_biomass_2020_j_index_v20240514').select('total_biomass')
zones_fc = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/zones').select('FIRST_zone', 'dsl')

# Convert zones to image and get metadata
zones_img = zones_fc.reduceToImage(['dsl'], ee.Reducer.min()).byte().updateMask(total_biomass.mask())
zone_codes = zones_fc.aggregate_array('dsl').add(4)
zone_dict = ee.Dictionary({'1': 'High Arctic', '2': 'Low Arctic', '3': 'Oro Arctic', '4': 'Pan Arctic'})

# Create area image
area = ee.Image.pixelArea()

# Combine zone and area images
zones_area = zones_img.addBands(area)

# =========================
# 2. ANALYSIS =============
# =========================

# 2.0 ----- CALCULATE AREA ACROSS BIOCLIMATE ZONES -----

zones_area_totals = []
for zone_code in zone_codes.getInfo():

    # Get zone name
    zone_name = zone_dict.get(ee.Number(zone_code).format('%.0f'))

    # Mask biomass to current code only
    current_zone = zones_area.updateMask(zones_img.eq(ee.Image.constant(zone_code)).selfMask())

    # Tidy mask
    current_zone = current_zone.mask(current_zone.mask().gt(0))

    # Reduce over masked biomass image
    current_zone_area = current_zone.reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=ee.Geometry.Polygon(params.export_region_coords, None, False),
        scale=1000,
        crs='EPSG:4326',
        maxPixels=1e13,
        tileScale=1)  # Default: 1

    # Convert to feature
    feat = ee.Feature(None, current_zone_area).set('dsl', zone_code).set('FIRST_zone', zone_name)

    zones_area_totals.append(feat)

# Export
task = ee.batch.Export.table.toDrive(
    collection=ee.FeatureCollection(zones_area_totals).select(['.*'], None, False),
    description='zone_area',
    folder='areas'
)
task.start()

raise Exception('stop')