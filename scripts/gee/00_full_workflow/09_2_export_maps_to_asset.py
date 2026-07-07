"""

DESCRIPTION: Export maps to asset

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

# Visualization options
palette = ['2a186c', '1b3d8f', '1566a5', '1596be', '30c8c2', '74e9bd', 'c5f9d0']  # cmocean.Haline[7]

# Reprojection options
crs = params.crs
scale = params.scale
biomass_palette = ['2a186c', '1b3d8f', '1566a5', '1596be', '30c8c2', '74e9bd', 'c5f9d0']  # cmocean.Haline[7]  # Biomass (g/m2)

# Model options
version = 'v20240514'
biomass_max = 3000
year = 2000
threshold_metric = 'j_index'

#  1.1 ----- READ IN DATA ----

total_continuous_mc = ee.ImageCollection('projects/arctic-biomass-mapping/assets/modeled_continuous/model_continuous_' + str(year) + '_total_' + version + '_mc')
woody_continuous_mc = ee.ImageCollection('projects/arctic-biomass-mapping/assets/modeled_continuous/model_continuous_' + str(year) + '_woody_' + version + '_mc')
total_probability_mc = ee.ImageCollection('projects/arctic-biomass-mapping/assets/modeled_binary/model_binary_' + str(year) + '_total_' + version + '_mc')
woody_probability_mc = ee.ImageCollection('projects/arctic-biomass-mapping/assets/modeled_binary/model_binary_' + str(year) + '_woody_' + version + '_mc')

if threshold_metric == 'j_index':
    total_probability_thresholds = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/model_ready_data/v20240514/mc/binary/thresholds/presence_thresholds_j_index_total_binary_v20240514')
    woody_probability_thresholds = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/model_ready_data/v20240514/mc/binary/thresholds/presence_thresholds_j_index_woody_binary_v20240514')
else:
    total_probability_thresholds = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/model_ready_data/v20240514/mc/binary/thresholds/presence_thresholds_f_meas_total_binary_v20240514')
    woody_probability_thresholds = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/model_ready_data/v20240514/mc/binary/thresholds/presence_thresholds_f_meas_woody_binary_v20240514')

# ============================
# 2. CALCULATE PERCENTILES ===
# ============================

#  2.0 ----- CONTINUOUS -----

total_continuous_percentiles = total_continuous_mc.reduce(ee.Reducer.percentile([2, 3, 50, 97, 98]))
woody_continuous_percentiles = woody_continuous_mc.reduce(ee.Reducer.percentile([2, 3, 50, 97, 98]))

#  2.1 ----- PROBABILITY -----

total_probability_percentiles = total_probability_mc.reduce(ee.Reducer.percentile([2, 3, 50, 97, 98]))
woody_probability_percentiles = woody_probability_mc.reduce(ee.Reducer.percentile([2, 3, 50, 97, 98]))

# =================================================
# 3. DERIVE PRESENCE/ABSENCE FROM PROBABILITIES ===
# =================================================

#  3.1 -----  APPLY PRESENCE THRESHOLDS -----

# Set up list of Monte Carlo iteration numbers as strings
mc_iterations = ee.List.sequence(1, 100, 1).map(lambda item: ee.Number(item).format('%.0f')).getInfo()

# Apply threshold - total
def apply_total_threshold(mc):
    probability = total_probability_mc.filter(ee.Filter.eq('system:index', mc)).first()
    threshold = total_probability_thresholds.filter(ee.Filter.eq('mc', ee.Number.parse(mc)))
    threshold = ee.Number(threshold.first().get('threshold')).multiply(100).round()
    presence_absence = probability.gte(threshold).rename('presence').byte()
    return presence_absence

total_binary = ee.ImageCollection([apply_total_threshold(mc) for mc in mc_iterations])

# Apply threshold - woody
def apply_woody_threshold(mc):
    probability = woody_probability_mc.filter(ee.Filter.eq('system:index', mc)).first()
    threshold = woody_probability_thresholds.filter(ee.Filter.eq('mc', ee.Number.parse(mc)))
    threshold = ee.Number(threshold.first().get('threshold')).multiply(100).round()
    presence_absence = probability.gte(threshold).rename('presence').byte()
    return presence_absence

woody_binary = ee.ImageCollection([apply_woody_threshold(mc) for mc in mc_iterations])

#  3.2 ----- COMPOSITE - MODE -----

total_presence_absence = total_binary.mode()
woody_presence_absence = woody_binary.mode()


# =========================
# 4. FINALIZE MAPS ========
# =========================

# 4.0 ----- REMOVE MASK OPACITY -----

# Continuous
total_continuous_percentiles = total_continuous_percentiles.updateMask(total_continuous_percentiles.mask().gt(0))
woody_continuous_percentiles = woody_continuous_percentiles.updateMask(woody_continuous_percentiles.mask().gt(0))

# Presence/absence
total_presence_absence = total_presence_absence.updateMask(total_presence_absence.mask().gt(0))
woody_presence_absence = woody_presence_absence.updateMask(woody_presence_absence.mask().gt(0))

# Probability
total_probability_percentiles = total_probability_percentiles.updateMask(total_probability_percentiles.mask().gt(0))
woody_probability_percentiles = woody_probability_percentiles.updateMask(woody_probability_percentiles.mask().gt(0))

# 4.1 ----- COMBINE -----

# Biomass
total_biomass_percentiles = total_continuous_percentiles.multiply(total_presence_absence).uint16()
woody_biomass_percentiles = woody_continuous_percentiles.multiply(total_presence_absence).multiply(woody_presence_absence).uint16()

# 4.2 ----- PARTITION -----

# Total biomass
total_biomass_median = total_biomass_percentiles.select('predicted_p50').rename('total_biomass')
total_biomass_lwr = ee.ImageCollection([total_biomass_percentiles.select('predicted_p2').rename('total_biomass_lwr'), total_biomass_percentiles.select('predicted_p3').rename('total_biomass_lwr')]).mean()
total_biomass_upr = ee.ImageCollection([total_biomass_percentiles.select('predicted_p97').rename('total_biomass_upr'), total_biomass_percentiles.select('predicted_p98').rename('total_biomass_upr')]).mean()

# Woody biomass
woody_biomass_median = woody_biomass_percentiles.select('predicted_p50').rename('woody_biomass')
woody_biomass_lwr = ee.ImageCollection([woody_biomass_percentiles.select('predicted_p2').rename('woody_biomass_lwr'), woody_biomass_percentiles.select('predicted_p3').rename('woody_biomass_lwr')]).mean()
woody_biomass_upr = ee.ImageCollection([woody_biomass_percentiles.select('predicted_p97').rename('woody_biomass_upr'), woody_biomass_percentiles.select('predicted_p98').rename('woody_biomass_upr')]).mean()

# Total probability
total_probability_lwr = ee.ImageCollection([total_probability_percentiles.select('prob_presence_p2').rename('total_probability_lwr'), total_probability_percentiles.select('prob_presence_p3').rename('total_probability_lwr')]).mean()
total_probability_upr = ee.ImageCollection([total_probability_percentiles.select('prob_presence_p97').rename('total_probability_upr'), total_probability_percentiles.select('prob_presence_p98').rename('total_probability_upr')]).mean()

# Woody probability
woody_probability_lwr = ee.ImageCollection([woody_probability_percentiles.select('prob_presence_p2').rename('woody_probability_lwr'), woody_probability_percentiles.select('prob_presence_p3').rename('woody_probability_lwr')]).mean()
woody_probability_upr = ee.ImageCollection([woody_probability_percentiles.select('prob_presence_p97').rename('woody_probability_upr'), woody_probability_percentiles.select('prob_presence_p98').rename('woody_probability_upr')]).mean()

# 4.3 ----- ROUND TOTAL UP TO MATCH WOODY -----
# This prevents instances where woody biomass > total biomass which is not ecologically sensible

total_biomass_median = total_biomass_median.max(woody_biomass_median)

# 4.4 ----- RE-COMBINE -----

# Biomass
total_biomass = total_biomass_median.addBands(total_biomass_lwr).addBands(total_biomass_upr)
woody_biomass = woody_biomass_median.addBands(woody_biomass_lwr).addBands(woody_biomass_upr)

# Probability
total_probability = total_probability_percentiles.select('prob_presence_p50').rename('total_probability').addBands(total_probability_lwr).addBands(total_probability_upr)
woody_probability = woody_probability_percentiles.select('prob_presence_p50').rename('woody_probability').addBands(woody_probability_lwr).addBands(woody_probability_upr)

# 4.5 ----- MAP -----

# Biomass

# Probability

# =========================
# 5. EXPORT ===============
# =========================

task = ee.batch.Export.image.toAsset(
    image=total_biomass.round().uint16(),  # Convert to int16 for easier storage
    description='total_biomass_' + str(year) + '_' + threshold_metric + '_' + version,
    assetId='projects/arctic-biomass-mapping/assets/modeled_final/total_biomass_' + str(year) + '_' + threshold_metric + '_' + version,
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    maxPixels=1e13
)
task.start()

task = ee.batch.Export.image.toAsset(
    image=woody_biomass.round().uint16(),  # Convert to int16 for easier storage
    description='woody_biomass_' + str(year) + '_' + threshold_metric + '_' + version,
    assetId='projects/arctic-biomass-mapping/assets/modeled_final/woody_biomass_' + str(year) + '_' + threshold_metric + '_' + version,
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    maxPixels=1e13
)
task.start()

task = ee.batch.Export.image.toAsset(
    image=total_probability.round().byte(),  # Convert to byte for easier storage
    description='total_probability_' + str(year) + '_' + threshold_metric + '_' + version,
    assetId='projects/arctic-biomass-mapping/assets/modeled_final/total_probability_' + str(year) + '_' + threshold_metric + '_' + version,
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    maxPixels=1e13
)
task.start()

task = ee.batch.Export.image.toAsset(
    image=woody_probability.round().byte(),  # Convert to byte for easier storage
    description='woody_probability_' + str(year) + '_' + threshold_metric + '_' + version,
    assetId='projects/arctic-biomass-mapping/assets/modeled_final/woody_probability_' + str(year) + '_' + threshold_metric + '_' + version,
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    maxPixels=1e13
)
task.start()

task = ee.batch.Export.image.toAsset(
    image=total_presence_absence.byte(),  # Convert to byte for easier storage
    description='total_presence_absence_' + str(year) + '_' + threshold_metric + '_' + version,
    assetId='projects/arctic-biomass-mapping/assets/modeled_final/total_presence_absence_' + str(year) + '_' + threshold_metric + '_' + version,
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    maxPixels=1e13
)
task.start()

task = ee.batch.Export.image.toAsset(
    image=woody_presence_absence.byte(),  # Convert to byte for easier storage
    description='woody_presence_absence_' + str(year) + '_' + threshold_metric + '_' + version,
    assetId='projects/arctic-biomass-mapping/assets/modeled_final/woody_presence_absence_' + str(year) + '_' + threshold_metric + '_' + version,
    region=ee.Geometry.Polygon(params.export_region_coords, None, False),  # Specify as rectangle to avoid gaps in output
    scale=scale,
    crs=crs,
    maxPixels=1e13
)
task.start()

raise Exception('stop')