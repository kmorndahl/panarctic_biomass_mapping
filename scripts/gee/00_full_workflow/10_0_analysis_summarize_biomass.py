"""

DESCRIPTION: Summarize biomass amount (Tg) and density (g/m2) across multiple subcategories:
bioclimate zones, CAVM vegetation community types; countries; Raynolds and Spawn project extents

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:
- Tile scales are adjusted for each dataset to facilitate exports
- For comparison to Raynolds and Spawn products:
  - We need to request data at 30m so that the mask is applied properly
  - Requesting at larger scales, the water/ice mask is applied at that larger scale and is inaccurate
  - In these instances, the disaggregation to 30 m starts at the beginning of the workflow, for example for the Raynolds product:
    - The original image is in density units, kg/m2
    - This image with ~8 km pixels is deconstructed into 30 m pixels, for example
      - 8 km pixel value = 0.7 kg/m2, each 30 m pixel that makes up this larger pixel will also be 0.7 kg/m2
      - So, this 8 km pixel would contain about 44,800,000 kg of biomass total (0.7 kg/m2 * 64,000,000 m2 pixel area)
        - Each 30 m pixel would contain about 630 kg (0.7 kg/m2 * 900 m2 pixel area)
        - Approximately 71,111 30 m pixels in an 11 km pixel
        - (630 kg / 30 m pixel) * 71,111 30 m pixels in an 8 km pixel = 44,799,930 (slightly off due to rounding)

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

# Visualization options
palette = ['2a186c', '1b3d8f', '1566a5', '1596be', '30c8c2', '74e9bd', 'c5f9d0']  # cmocean.Haline[7]

# Reprojection options
crs = params.crs
canopy_height_threshold = params.canopy_height_threshold
scale = params.scale  # Scale for reprojecting predictors
biomass_palette = ['2a186c', '1b3d8f', '1566a5', '1596be', '30c8c2', '74e9bd', 'c5f9d0']  # Biomass (g/m2)

# Model options
version = 'v20240514'
year = 2000

# Mask options
mask_comparison = 'masked'  # Set to 'masked' to mask Orndahl/Spawn/Raynolds to exact same extent, set to 'unmasked' to use original resolution and extent of each
output_folder = 'biomass_summaries'

# 1.1 ----- READ IN DATA -----

total_biomass = ee.Image('projects/arctic-biomass-mapping/assets/modeled_final/total_biomass_' + str(year) + '_j_index_v20240514')
woody_biomass = ee.Image('projects/arctic-biomass-mapping/assets/modeled_final/woody_biomass_' + str(year) + '_j_index_v20240514')
canopy_height = ee.ImageCollection("projects/meta-forest-monitoring-okw37/assets/CanopyHeight").mosaic()
countries = ee.FeatureCollection("WM/geoLab/geoBoundaries/600/ADM0").filter(ee.Filter.inList('shapeGroup', ['USA', 'CAN', 'RUS', 'GRL', 'ISL', 'NOR', 'SWE', 'FIN']))
zones_fc = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/zones').select('FIRST_zone', 'dsl')
cavm = ee.Image('projects/arctic-biomass-mapping/assets/ROIs/cavm_raster')
cavm_legend = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/cavm_raster_legend')
total_biomass_raynolds = ee.Image('projects/arctic-biomass-mapping/assets/compare/aga_circumpolar_avhrr_biomass_2010_epsg3571')
total_biomass_spawn = ee.ImageCollection("NASA/ORNL/biomass_carbon_density/v1").first().select(['agb', 'bgb'])

# Create reducer
# Use weighted reducers so fractional pixels are represented properly
# sum = total biomass per zone/region (Tg)
# mean = biomass density per zone/region (g/m2)
reducer = ee.Reducer.mean().combine(
    reducer2=ee.Reducer.sum(),
    sharedInputs=True
)

# =========================
# 2. TIDY DATA ============
# =========================

# 2.0 ----- CALCULATE TOTAL GRAMS PER PIXEL -----

# Area image
# Area will be calculated at whatever scale is specified for reduceRegion/export
area = ee.Image.pixelArea()

# Multiply by pixel area to get total grams in each pixel
total_biomass_g = total_biomass.multiply(area).regexpRename('$', '_g')
woody_biomass_g = woody_biomass.multiply(area).regexpRename('$', '_g')

# 2.1 ----- WOODY PLANT DOMINANCE -----

# Calculate woody percent
woody_percent_mask = woody_biomass.select('woody_biomass').add(total_biomass.select('total_biomass')).eq(0).Not()  # Mask if both total and woody biomass are 0
woody_percent_vegetated = (woody_biomass.select('woody_biomass')
                                        .divide(total_biomass.select('total_biomass'))
                                        .multiply(100)
                                        .updateMask(woody_percent_mask)
                                        .select('woody_biomass')
                                        .rename('woody_percent_vegetated'))  # Original - mask out where plant and woody plant biomass are both zero
woody_percent = (woody_biomass.select('woody_biomass')
                              .divide(total_biomass.select('total_biomass'))
                              .multiply(100)
                              .select('woody_biomass')
                              .rename('woody_percent'))  # New option 'set0' - don't mask, where plant and woody plant biomass are both zero, woody percent is zero

# 2.2 ----- MASK -----

# Mask by canopy height
canopy_height = canopy_height.unmask(0)  # Fill currently masked areas with low value to ensure we include Greeland, barrier islands etc.
canopy_height_mask = canopy_height.lt(canopy_height_threshold).selfMask()
total_biomass_g = total_biomass_g.updateMask(canopy_height_mask)
woody_biomass_g = woody_biomass_g.updateMask(canopy_height_mask)

# 2.3 ----- COMBINE -----

biomass_g = total_biomass_g.addBands(woody_biomass_g).addBands(woody_percent).addBands(woody_percent_vegetated)

# =========================
# 3. SUMMARIZE ============
# =========================

# 3.0 ----- COUNTRIES -----

# Summarize over countries
biomass_g_countries = biomass_g.clip(countries).reduceRegions(
    collection=countries,
    reducer=reducer,
    scale=scale,
    crs=crs,
    tileScale=1
)

print('Biomass summarized over countries', biomass_g_countries.first())

# Export
task = ee.batch.Export.table.toDrive(
    collection=biomass_g_countries.select(['.*'], None, False),
    description='biomass_country_summary_mask' + str(canopy_height_threshold) + '_' + str(year) + '_' + version,
    folder=output_folder
)
task.start()

# Union countries to create Pan Arctic
pan_arctic = (countries.union()
                       .set('shapeGroup', 'PAN')
                       .set('shapeName', 'Pan Arctic')
                       .set('shapeType', 'ADM0'))

# Tidy mask
biomass_g_pan = biomass_g.mask(biomass_g.mask().gt(0))

# Reduce over Pan Arctic
# Must be done separately to avoid out of memory errors
pan_arctic_summary = biomass_g_pan.clip(pan_arctic).reduceRegions(
    collection=pan_arctic,
    reducer=reducer,
    scale=scale,
    crs=crs,
    tileScale=8
)

# Export
task = ee.batch.Export.table.toDrive(
    collection=pan_arctic_summary.select(['.*'], None, False),
    description='biomass_country_summary_PAN_mask' + str(canopy_height_threshold) + '_' + str(year) + '_' + version,
    folder=output_folder
)
task.start()

# 3.1 ----- BIOCLIMATE ZONES -----

# Tidy zones
zones_img = zones_fc.reduceToImage(['dsl'], ee.Reducer.min()).byte()

# Create arctic image
arctic_img = zones_img.selfMask()

# Get zone information
zone_codes = zones_fc.aggregate_array('dsl').add(4)
zone_dict = ee.Dictionary({'1': 'High Arctic', '2': 'Low Arctic', '3': 'Oro Arctic', '4': 'Pan Arctic'})

# Summarize over bioclimate zones
biomass_g_zones = []
for zone_code in zone_codes.getInfo():

    # Get zone name
    zone_name = zone_dict.get(ee.Number(zone_code).format('%.0f'))

    # Mask biomass to current code only
    biomass_g_zone_code = ee.Image(ee.Algorithms.If(
        ee.Number(zone_code).eq(4),
        biomass_g.updateMask(arctic_img),
        biomass_g.updateMask(zones_img.eq(ee.Image.constant(zone_code)).selfMask())))

    # Tidy mask
    biomass_g_zone_code = biomass_g_zone_code.mask(biomass_g_zone_code.mask().gt(0))

    # Reduce over masked biomass image
    biomass_g_zone_code_reduced = biomass_g_zone_code.reduceRegion(
        reducer=reducer,
        geometry=ee.Geometry.Polygon(params.export_region_coords, None, False),
        scale=scale,
        crs=crs,
        maxPixels=1e13,
        tileScale=0.5)  # Default: 1

    # Convert to feature
    feat = ee.Feature(None, biomass_g_zone_code_reduced).set('dsl', zone_code).set('FIRST_zone', zone_name)

    biomass_g_zones.append(feat)

print('Biomass summarized over bioclimate zones', ee.FeatureCollection(biomass_g_zones).first())

# Export
task = ee.batch.Export.table.toDrive(
    collection=ee.FeatureCollection(biomass_g_zones).select(['.*'], None, False),
    description='biomass_zone_summary_mask' + str(canopy_height_threshold) + '_' + str(year) + '_' + version,
    folder=output_folder
)
task.start()

# 3.2 ----- CAVM VEGETATION COMMUNITY TYPES (FINE) -----

# Tidy data
cavm_codes_fine = cavm_legend.aggregate_array('raster_code')

# Pre-fetch CAVM legend as a Python dict before the loop.
# Each feature in cavm_legend has 'raster_code' and 'short_description' properties; fetching
# all at once avoids a separate GEE filter expression per code inside the loop.
cavm_dict = {
    feat['properties']['raster_code']: feat['properties']['short_description']
    for feat in cavm_legend.getInfo()['features']
}

# Summarize over CAVM vegetation community types (fine)
biomass_g_cavm_fine = []
for code in cavm_codes_fine.getInfo():

    # Get vegetation description
    veg_desc = cavm_dict[code]

    # Restrict CAVM image to current code only
    cavm_code = cavm.eq(ee.Image.constant(code)).selfMask()

    # Mask biomass to current code only
    biomass_g_cavm_code = biomass_g.updateMask(cavm_code)

    # Tidy mask
    biomass_g_cavm_code = biomass_g_cavm_code.mask(biomass_g_cavm_code.mask().gt(0))

    # Reduce over masked biomass image
    biomass_g_cavm_code_reduced = biomass_g_cavm_code.reduceRegion(
        reducer=reducer,
        geometry=ee.Geometry.Polygon(params.export_region_coords, None, False),
        scale=scale,
        crs=crs,
        maxPixels=1e13,
        tileScale=0.2)  # Default: 0.5

    # Convert to feature
    feat = ee.Feature(None, biomass_g_cavm_code_reduced).set('veg_code', code).set('veg_desc', veg_desc)

    biomass_g_cavm_fine.append(feat)

print('Biomass summarized over CAVM vegetation community types (fine)', ee.FeatureCollection(biomass_g_cavm_fine).first())

# Export
task = ee.batch.Export.table.toDrive(
    collection=ee.FeatureCollection(biomass_g_cavm_fine).select(['.*'], None, False),
    description='biomass_cavm_fine_summary_mask' + str(canopy_height_threshold) + '_' + str(year) + '_' + version,
    folder=output_folder
)
task.start()

# 3.3 ----- CAVM VEGETATION COMMUNITY TYPES (COARSE) -----

# Convert CAVM codes to coarse scale
cavm_coarse = cavm.divide(10).floor()  # Convert to coarser aggregation level

def to_coarse(item):
    return ee.Number.parse(item).divide(10).floor()

cavm_codes_coarse = cavm_codes_fine.map(to_coarse).distinct()  # Convert to coarser aggregation level

# Summarize over CAVM vegetation community types (coarse)
biomass_g_cavm_coarse = []
for code in cavm_codes_coarse.getInfo():

    # Restrict CAVM image to current code only
    cavm_code = cavm_coarse.eq(ee.Image.constant(code)).selfMask()

    # Mask biomass to current code only
    biomass_g_cavm_code = biomass_g.updateMask(cavm_code)

    # Tidy mask
    biomass_g_cavm_code = biomass_g_cavm_code.mask(biomass_g_cavm_code.mask().gt(0))

    # Reduce over masked biomass image
    biomass_g_cavm_code_reduced = biomass_g_cavm_code.reduceRegion(
        reducer=reducer,
        geometry=ee.Geometry.Polygon(params.export_region_coords, None, False),
        scale=scale,
        crs=crs,
        maxPixels=1e13,
        tileScale=0.5)  # Default: 1

    # Convert to feature
    feat = ee.Feature(None, biomass_g_cavm_code_reduced).set('veg_code', code)

    biomass_g_cavm_coarse.append(feat)

print('Biomass summarized over CAVM vegetation community types (coarse)', ee.FeatureCollection(biomass_g_cavm_coarse).first())

# Export
task = ee.batch.Export.table.toDrive(
    collection=ee.FeatureCollection(biomass_g_cavm_coarse).select(['.*'], None, False),
    description='biomass_cavm_coarse_summary_mask' + str(canopy_height_threshold) + '_' + str(year) + '_' + version,
    folder=output_folder
)
task.start()

# =========================
# 4. COMPARISONS ==========
# =========================

# 4.0 ----- SET UP ROI -----

# Restrict to High and Low Arctic only to match extent of Raynolds map
zones_img_hl = zones_img.mask(zones_img.lt(3))
arctic_img_hl = zones_img_hl.selfMask()
zone_codes_hl = ee.List([1, 2, 4])

# 4.0 ----- MASK -----

# Mask permanent ice and water from Raynolds and Spawn maps to match masking in our maps
mask_raynolds_total = total_biomass_g.select('total_biomass_g').mask()
mask_spawn_total = total_biomass_g.select('total_biomass_g').mask()

# 4.1 ----- RAYNOLDS -----

# Rescale Raynolds data
total_biomass_raynolds = total_biomass_raynolds.multiply(1000)  # Convert to g/m2
total_biomass_raynolds_g = total_biomass_raynolds.multiply(area).rename('total_biomass_g')  # Multiply by pixel area to get total grams in each pixel

# Mask each image to the other so that extents match
total_biomass_g_raynolds_extent = total_biomass_g.updateMask(mask_raynolds_total)
total_biomass_raynolds_g = total_biomass_raynolds_g.updateMask(mask_raynolds_total)

# Create unique band names
biomass_compare_raynolds = total_biomass_g_raynolds_extent.regexpRename('$', '_raynolds_extent').addBands(total_biomass_raynolds_g.regexpRename('$', '_raynolds'))

# Summarize over bioclimate zones
biomass_raynolds_summary = []
for zone_code in zone_codes_hl.getInfo():

    # Get zone name
    zone_name = zone_dict.get(ee.Number(zone_code).format('%.0f'))

    # Mask biomass to current code only
    biomass_g_raynolds_code = ee.Image(ee.Algorithms.If(
        ee.Number(zone_code).eq(4),
        biomass_compare_raynolds.updateMask(arctic_img_hl),
        biomass_compare_raynolds.updateMask(zones_img_hl.eq(ee.Image.constant(zone_code)).selfMask())))

    # Tidy mask
    biomass_g_raynolds_code = biomass_g_raynolds_code.mask(biomass_g_raynolds_code.mask().gt(0))

    # Reduce over masked biomass image
    biomass_g_raynolds_code_reduced = biomass_g_raynolds_code.reduceRegion(
        reducer=reducer,
        geometry=ee.Geometry.Polygon(params.export_region_coords, None, False),
        scale=scale,
        crs=crs,
        maxPixels=1e13,
        tileScale=0.2)  # Default: 1

    # Convert to feature
    feat = ee.Feature(None, biomass_g_raynolds_code_reduced).set('dsl', zone_code).set('FIRST_zone', zone_name)

    biomass_raynolds_summary.append(feat)

print('Biomass summarized over bioclimate zones and compared to Raynolds et al.', ee.FeatureCollection(biomass_raynolds_summary).first())

# Export
task = ee.batch.Export.table.toDrive(
    collection=ee.FeatureCollection(biomass_raynolds_summary).select(['.*'], None, False),
    description='biomass_compare_raynolds_' + mask_comparison + '_summary_mask' + str(canopy_height_threshold) + '_' + str(year) + '_' + version,
    folder=output_folder
)
task.start()

# 4.2 ----- SPAWN -----

# Rescale data
total_biomass_spawn = total_biomass_spawn.divide(0.492)  # Convert to AGB from AGBC
total_biomass_spawn = total_biomass_spawn.multiply(100)  # Convert from Mg/ha to g/m2

# Multiply by pixel area to get total grams in each pixel
total_biomass_spawn_g = total_biomass_spawn.multiply(area).rename(['agb_g', 'bgb_g'])

# Mask each image to the other so that extents match
total_biomass_g_spawn_extent = total_biomass_g.updateMask(mask_spawn_total)
total_biomass_spawn_g = total_biomass_spawn_g.updateMask(mask_spawn_total)

# Create unique band names
biomass_compare_spawn = (total_biomass_g_spawn_extent.regexpRename('$', '_spawn_extent')
                                                     .addBands(total_biomass_spawn_g.regexpRename('$', '_spawn')))

# Summarize over bioclimate zones
biomass_spawn_summary = []
for zone_code in zone_codes_hl.getInfo():

    # Get zone name
    zone_name = zone_dict.get(ee.Number(zone_code).format('%.0f'))

    # Mask biomass to current code only
    biomass_g_spawn_code = ee.Image(ee.Algorithms.If(
        ee.Number(zone_code).eq(4),
        biomass_compare_spawn.updateMask(arctic_img_hl),
        biomass_compare_spawn.updateMask(zones_img_hl.eq(ee.Image.constant(zone_code)).selfMask())))

    # Tidy mask
    biomass_g_spawn_code = biomass_g_spawn_code.mask(biomass_g_spawn_code.mask().gt(0))

    # Reduce over masked biomass image
    biomass_g_spawn_code_reduced = biomass_g_spawn_code.reduceRegion(
        reducer=reducer,
        geometry=ee.Geometry.Polygon(params.export_region_coords, None, False),
        scale=scale,
        crs=crs,
        maxPixels=1e13,
        tileScale=0.5)  # Default: 1

    # Convert to feature
    feat = ee.Feature(None, biomass_g_spawn_code_reduced).set('dsl', zone_code).set('FIRST_zone', zone_name)

    biomass_spawn_summary.append(feat)

print('Biomass summarized over bioclimate zones and compared to Spawn et al.', ee.FeatureCollection(biomass_spawn_summary).first())

# Export
task = ee.batch.Export.table.toDrive(
    collection=ee.FeatureCollection(biomass_spawn_summary).select(['.*'], None, False),
    description='biomass_compare_spawn_' + mask_comparison + '_summary_mask' + str(canopy_height_threshold) + '_' + str(year) + '_' + version,
    folder=output_folder
)
task.start()

raise Exception('stop')