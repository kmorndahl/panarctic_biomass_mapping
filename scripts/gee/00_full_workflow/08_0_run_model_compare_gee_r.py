"""

DESCRIPTION: Confirm that models applied in GEE exactly match results from R

AUTHOR: Kathleen Orndahl
DATE: 11-10-2024

NOTES:
  - For binary, R models are probability forests
  - Therefore, need to format these as continuous trees with probabilities when converting to GEE tree format
  - https://stackoverflow.com/questions/62806074/how-to-get-the-same-prediction-probability-and-class-in-a-random-forest

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
response_type = 'binary'
ds_type = 'total'
mc = 1

# Cal/val extraction options
calval_buffer_diam = 90  # 90m x 90m window to extract site level data
calval_buffer = ee.Number(calval_buffer_diam).divide(2)
band_name = ee.Algorithms.If(
    ee.String(response_type).equals('binary'),
    'raw_probability',
    'predicted'
)
export_bands = ee.Algorithms.If(
    ee.String(response_type).equals('binary'),
    ['site_code', 'prob_presence', 'prob_absence'],
    ['site_code', 'predicted']
)

# 1.1 ----- READ IN DATA -----

calval_data = ee.FeatureCollection(
    "projects/arctic-biomass-mapping/assets/model_ready_data/" + version + "/mc/test_gee_r_mod_match/final_data_" + ds_type + "_" + response_type + "_" + version + "_" + str(mc)
)  # Full dataset

# print('Training data: ', calval_data.getInfo())

# 1.2 ----- FUNCTIONS -----

def getProbPresence(feat):
    prob_list = ee.Array(feat.get('raw_probability')).toList()
    prob_presence = ee.Number(prob_list.reduce(ee.Reducer.mean().unweighted()))
    return feat.set('prob_presence', prob_presence).set('prob_absence', ee.Number(1).subtract(prob_presence))

# =========================
# 2. MODEL ================
# =========================

# Get model string from Google Storage
mod_str = ee.List([
    ee.Blob(
        'gs://arctic_biomass_mapping_models/' + version + '/mc/' + response_type + '/final_trees/' + ds_type + '_' + response_type + '_formatted_gee_short_' + version + '_' + str(mc) + '.txt'
    ).string()
])
# print('Model string:', mod_str.getInfo())
# print('Model string length: ', mod_str.join("\n").length().getInfo())

# Convert to GEE classifier
mod = ee.Algorithms.If(
    ee.String(response_type).equals('binary'),
    ee.Classifier.decisionTreeEnsemble(mod_str).setOutputMode('RAW'),
    ee.Classifier.decisionTreeEnsemble(mod_str)
)
mod = ee.Classifier(mod)
# print('Model:', mod.getInfo())

# Get information about the trained model
# print('Results of trained model', mod.explain().getInfo())

# Predict on training data to compare to R results
predictions = calval_data.classify(mod, band_name)

# If binary, further tidy predictions
# https://gis.stackexchange.com/questions/447452/predict-classification-probability-of-pretrained-random-forest-in-gee
predictions = ee.Algorithms.If(
    ee.String(response_type).equals('binary'),
    predictions.map(getProbPresence),
    predictions
)
predictions = ee.FeatureCollection(predictions)

# =========================
# 3. EXPORT ===============
# =========================

# Select predictors for export
predictions = predictions.select(export_bands)
# print('Final predictions: ', predictions.getInfo())

# Export
task = ee.batch.Export.table.toDrive(
    collection=predictions,
    description='GEE_predictions_' + response_type + '_rf_' + ds_type + '_' + version + '_' + str(mc),
    folder='compare_gee_r_model_predictions',
    fileFormat='CSV'
)
task.start()

raise Exception('stop')