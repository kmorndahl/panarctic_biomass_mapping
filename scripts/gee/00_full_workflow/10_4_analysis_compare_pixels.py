"""

DESCRIPTION: Sample pixels from 30m, 300m and 8km products for comparison

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

# Sampling options
prop = 1  # Proportion of pixels
perc = prop * 100

# Reprojection options
crs = params.crs

# Output options
output_folder = 'compare_pixels'

# 1.1 ----- READ IN DATA -----

zones_fc = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/zones').select('FIRST_zone', 'dsl')
total_biomass = ee.Image("projects/arctic-biomass-mapping/assets/modeled_final/total_biomass_2020_j_index_v20240514")
total_biomass_spawn = ee.ImageCollection("NASA/ORNL/biomass_carbon_density/v1").first().select('agb')
total_biomass_raynolds = ee.Image('projects/arctic-biomass-mapping/assets/compare/aga_circumpolar_avhrr_biomass_2010_epsg3571')

# 1.2 ----- TIDY DATA -----

# Restrict to high and low Arctic
zones = zones_fc.filter(ee.Filter.neq('FIRST_zone', 'Oro Arctic')).union()

# Rescale data to units of g/m2
total_biomass = total_biomass.select('total_biomass')  # Get median value only
total_biomass_spawn = total_biomass_spawn.divide(0.492)  # Convert to AGB from AGBC
total_biomass_spawn = total_biomass_spawn.multiply(100)  # Convert from Mg/ha to g/m2
total_biomass_raynolds = total_biomass_raynolds.multiply(1000)  # Convert to g/m2

# Mask permanent ice and water from Raynolds and Spawn maps to match masking in our maps
mask = total_biomass.mask()
total_biomass_raynolds = total_biomass_raynolds.updateMask(mask)
total_biomass_spawn = total_biomass_spawn.updateMask(mask)

# Get metadata
spawn_count = total_biomass_spawn.reduceRegion(reducer=ee.Reducer.count(), geometry=zones.geometry(), maxPixels=1e13)
raynolds_count = total_biomass_raynolds.reduceRegion(reducer=ee.Reducer.count())
print('Total number of pixels, Spawn: ', spawn_count)
print('Total number of pixels, Raynolds: ', raynolds_count)
print('Orndahl pixel scale: ', total_biomass.projection().nominalScale())
print('Spawn pixel scale: ', total_biomass_spawn.projection().nominalScale())
print('Raynolds pixel scale: ', total_biomass_raynolds.projection().nominalScale())

# =========================
# 2. ANALYSIS =============
# =========================

# 2.0 ----- COMBINE DATA -----

biomass_gm2_spawn = total_biomass.rename('biomass_density_gm2').addBands(total_biomass_spawn.rename('biomass_density_gm2_spawn'))
biomass_gm2_raynolds = total_biomass.rename('biomass_density_gm2').addBands(total_biomass_raynolds.rename('biomass_density_gm2_raynolds'))

# 2.1 ----- SAMPLE DATA -----

samples_gm2_spawn = biomass_gm2_spawn.sample(
    region=zones.geometry(),
    scale=total_biomass_spawn.projection().nominalScale(),
    projection=crs,
    dropNulls=True,
    geometries=False,
    factor=prop,
    seed=65,
    tileScale=1,
)

samples_gm2_raynolds = biomass_gm2_raynolds.sample(
    region=zones.geometry(),
    scale=total_biomass_raynolds.projection().nominalScale(),
    projection=crs,
    dropNulls=True,
    geometries=False,
    factor=prop,
    seed=65,
    tileScale=1
)

# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.table.toDrive(
    collection=samples_gm2_spawn,
    description='spawn_compare_pixels_avg_gm2_perc' + str(perc),
    folder=output_folder
)
task.start()

task = ee.batch.Export.table.toDrive(
    collection=samples_gm2_raynolds,
    description='raynolds_compare_pixels_avg_gm2_perc' + str(perc),
    folder=output_folder
)
task.start()

raise Exception('stop')