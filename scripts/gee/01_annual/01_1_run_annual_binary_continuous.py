"""

DESCRIPTION: Export binary, continuous, and probability predictions — MC-outer variant

AUTHOR: Kathleen Orndahl
DATE: 2026-07-06

NOTES:

- Variant of 01_1_run_annual_binary_continuous_agg_years.py that uses the MC-outer loop
  pattern for all three output types.

- ROOT CAUSE OF OOM (continuous):
  In the year-outer script, run_continuous_mc_ic_map() is called once per year, creating
  100 blob references per year. With 40 years stacked, the export graph contains
  100 × 40 = 4,000 blob references. GEE inlines blob string content per tile evaluation
  (~2.55 MB/blob × 4,000 = ~10 GB), far exceeding the per-tile memory limit.
  Binary exports succeed because blob files are much smaller (~42 MB total at 100 MC).

- MC-OUTER FIX:
  run_continuous_mc() and run_binary_mc() loop over MC in Python first.
  Within each MC iteration, a single ee.Blob / ee.Classifier Python object is created and
  reused across all year classify() calls. GEE's graph deduplicates shared Python objects,
  so the export graph contains exactly 100 blob references regardless of year count.

- SHARED IC FOR BINARY AND PROBABILITY:
  run_binary_mc() returns raw prob_presence × 100. The same IC is passed to both
  aggregate_binary_mc() (thresholding + mode) and aggregate_probability_mc()
  (percentile summary) so the binary model runs only once per ds_type.

- Output band names:
    Continuous:  {ds}_biomass_{year}, {ds}_biomass_lwr_{year}, {ds}_biomass_upr_{year}
    Binary:      {ds}_presence_absence_{year}
    Probability: {ds}_probability_{year}, {ds}_probability_lwr_{year}, {ds}_probability_upr_{year}

LOOP STRUCTURE:

  feature (outer) × pre-build predictors_by_year dict × ds_type × {cont MC IC, binary MC IC}

"""

# =========================
# 0. PATH SETUP ===========
# =========================

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import ee
from utils import params

try:
    ee.Initialize(project=params.ee_project_production)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=params.ee_project_production)

# =========================
# 1. SET UP ===============
# =========================

# 1.0 ----- PARAMETERS -----

# CCDC options
ccdc_version = 'v20240207'
ccdc_start_year = 1984
ccdc_end_year = 2023

# CCDC processing options
extrapolate_max_days = 120
segment_find_strategy = 'next'

# Model options
model_version = 'v20240514'
threshold_metric = 'j_index'  # 'j_index' or 'f_meas'
tree_dir = 'final_trees_shortened'

# Monte Carlo options
mc_list = list(range(1, 101))

# Auxiliary predictor options
texture_bands = ['spectral_NDVI_peakSummer']
texture_metrics = ['.*_var']
texture_radius = 1

# Reprojection options
crs = params.crs
scale = params.scale

# Input: agg-years Phase 1 reflectance collection (one image per feature, years as bands)
refl_version = '20260705'
refl_dir = 'modeled_reflectance'

# Output options
out_version = '20260705'
output_dir = 'annual_comparison'
output_base_path = 'projects/arctic-biomass-change/assets/' + output_dir + '/'
continuous_dir = 'continuous'
binary_dir = 'binary'
probability_dir = 'probability'

# Product types
ds_types = ['total', 'woody']

# Probability export (OPTIONAL)
export_probability = False

# Subset features (OPTIONAL)
subset_features = None  # None = process all features

# 1.1 ----- MODULES -----

from utils import pipeline
from utils import misc as fun_misc

# 1.2 ----- YEAR-INDEPENDENT DATA -----

print('Loading static predictors and CCDC data...')
static = pipeline.load_static_predictors()
ccdc_data = pipeline.load_ccdc_data(ccdc_version, ccdc_start_year, ccdc_end_year)

# Probability thresholds keyed by ds_type
if threshold_metric == 'j_index':
    probability_thresholds = {
        'total': ee.FeatureCollection(
            'projects/arctic-biomass-mapping/assets/model_ready_data/v20240514/mc/binary/thresholds/'
            'presence_thresholds_j_index_total_binary_v20240514'
        ),
        'woody': ee.FeatureCollection(
            'projects/arctic-biomass-mapping/assets/model_ready_data/v20240514/mc/binary/thresholds/'
            'presence_thresholds_j_index_woody_binary_v20240514'
        )
    }
else:
    probability_thresholds = {
        'total': ee.FeatureCollection(
            'projects/arctic-biomass-mapping/assets/model_ready_data/v20240514/mc/binary/thresholds/'
            'presence_thresholds_f_meas_total_binary_v20240514'
        ),
        'woody': ee.FeatureCollection(
            'projects/arctic-biomass-mapping/assets/model_ready_data/v20240514/mc/binary/thresholds/'
            'presence_thresholds_f_meas_woody_binary_v20240514'
        )
    }

# 1.3 ----- CREATE OUTPUT COLLECTIONS -----

fun_misc.createDirectoryIfNotExists(output_base_path + continuous_dir)
fun_misc.createDirectoryIfNotExists(output_base_path + binary_dir)

for ds in ds_types:
    fun_misc.createCollectionIfNotExists(
        output_base_path + continuous_dir + '/' + ds + '_continuous_' + out_version
    )
    fun_misc.createCollectionIfNotExists(
        output_base_path + binary_dir + '/' + ds + '_binary_' + out_version
    )

if export_probability:
    fun_misc.createDirectoryIfNotExists(output_base_path + probability_dir)
    for ds in ds_types:
        fun_misc.createCollectionIfNotExists(
            output_base_path + probability_dir + '/' + ds + '_probability_' + out_version
        )

# 1.4 ----- LOAD PHASE 1 REFLECTANCE COLLECTION -----

# Use aggregate_array() instead of toList().getInfo() — each image has ~1,400 bands, so
# serializing full image objects for 1,000 features exceeds the GEE user memory limit.
refl_ic = ee.ImageCollection(output_base_path + refl_dir + '/' + refl_dir + '_' + refl_version)
img_ids     = refl_ic.aggregate_array('system:id').getInfo()
feature_ids = refl_ic.aggregate_array('feature_id').getInfo()
year_starts = refl_ic.aggregate_array('year_start').getInfo()
year_ends   = refl_ic.aggregate_array('year_end').getInfo()

images_by_feature = {
    fid: {
        'id': img_id,
        'properties': {'feature_id': fid, 'year_start': ys, 'year_end': ye}
    }
    for fid, img_id, ys, ye in zip(feature_ids, img_ids, year_starts, year_ends)
}

feature_id_list = sorted(images_by_feature.keys())

if subset_features is not None:
    feature_id_list = feature_id_list[:subset_features]

# 1.5 ----- REPORT -----

print(f'Total features available: {len(images_by_feature)}')
print(f'Features to process (this run): {len(feature_id_list)}')

# =========================
# 2. PREDICTOR SHORTENING =
# =========================

_band_rename_cache = None  # (old_names, new_names); built on first call

def rename_predictors_short(predictors):
    """Rename predictor bands to shortened forms matching final_trees_shortened."""
    global _band_rename_cache
    if _band_rename_cache is None:
        old_names = predictors.bandNames().getInfo()
        new_names = [fun_misc.shortenPredName(b) for b in old_names]
        _band_rename_cache = (old_names, new_names)
    return predictors.rename(_band_rename_cache[1])

# ==============================
# 3. OUTER LOOP — FEATURES =====
# ==============================

for feature_id in feature_id_list:

    print(f'===== FEATURE: {feature_id} =====')

    feat_img_info = images_by_feature[feature_id]
    feat_img = ee.Image(feat_img_info['id'])

    year_start     = int(feat_img_info['properties']['year_start'])
    year_end       = int(feat_img_info['properties']['year_end'])
    year_list_feat = list(range(year_start, year_end + 1))

    region = feat_img.geometry()

    # -------------------------------------------------------------------------
    # 3.1 PRE-BUILD PREDICTOR IMAGES FOR ALL YEARS
    # Build all year-specific predictor images BEFORE the MC loop so that each
    # predictor is a single Python object. When the MC outer loop references
    # predictors_by_year[year] 100 times, GEE deduplicates the shared object and
    # serializes each predictor graph once — not once per MC iteration.
    # -------------------------------------------------------------------------

    print(f'  Pre-building predictors for {len(year_list_feat)} years...')
    predictors_by_year = {}

    for year in year_list_feat:

        modeled_reflectance = (feat_img
                               .select('.*_' + str(year))
                               .regexpRename('_' + str(year) + '$', ''))

        tree_cover = pipeline.get_tree_cover(year)
        aux = pipeline.compute_aux_predictors(
            modeled_reflectance,
            ccdc_data['segments'],
            year,
            crs,
            scale,
            texture_bands,
            texture_metrics,
            texture_radius,
            segment_find_strategy,
            extrapolate_max_days
        )
        spectral = pipeline.compute_spectral_predictors(
            modeled_reflectance, ccdc_data['seasonal_doy']
        )
        predictors = pipeline.assemble_predictors(
            spectral, static, tree_cover, aux['texture'], crs, scale
        )
        predictors_by_year[year] = rename_predictors_short(predictors)

    # -------------------------------------------------------------------------
    # 3.2 MODEL EXPORTS — MC-OUTER FOR ALL OUTPUT TYPES
    #
    # One ee.Blob / ee.Classifier per MC is shared across all year classify() calls,
    # keeping the export graph at exactly 100 blob references regardless of year count.
    #
    # Binary and probability share one IC (run_binary_mc called once per ds_type);
    # aggregate_binary_mc and aggregate_probability_mc both consume it.
    # -------------------------------------------------------------------------

    years_processed = year_list_feat

    for ds in ds_types:

        print(f'  Building MC graphs for {ds} ({len(year_list_feat)} years)...')

        # ----- Continuous (MC-outer) -----

        cont_mc_ic = pipeline.run_continuous_mc(
            ds, model_version, mc_list, predictors_by_year, year_list_feat, tree_dir
        )
        cont_stack = pipeline.aggregate_continuous_mc(cont_mc_ic, ds, year_list_feat)
        cont_stack = (cont_stack
                      .clip(region)
                      .set('feature_id', feature_id,
                           'years_processed', years_processed,
                           'n_years', len(years_processed))
                      .updateMask(cont_stack.mask().gt(0)))

        cont_asset_path = (output_base_path + continuous_dir + '/'
                           + ds + '_continuous_' + out_version + '/'
                           + str(feature_id))
        ee.batch.Export.image.toAsset(
            image=cont_stack,
            description=ds + '_continuous_' + str(feature_id),
            assetId=cont_asset_path,
            region=region,
            scale=scale,
            crs=crs,
            maxPixels=1e13
        ).start()
        print(f'  Submitted continuous: {cont_asset_path}')

        # ----- Binary and probability — shared MC IC (MC-outer) -----

        binary_mc_ic = pipeline.run_binary_mc(
            ds, model_version, mc_list, predictors_by_year, year_list_feat, tree_dir
        )

        bin_stack = pipeline.aggregate_binary_mc(binary_mc_ic, probability_thresholds[ds])
        bin_stack = (bin_stack
                     .clip(region)
                     .set('feature_id', feature_id,
                          'years_processed', years_processed,
                          'n_years', len(years_processed))
                     .updateMask(bin_stack.mask().gt(0)))

        bin_asset_path = (output_base_path + binary_dir + '/'
                          + ds + '_binary_' + out_version + '/'
                          + str(feature_id))
        ee.batch.Export.image.toAsset(
            image=bin_stack,
            description=ds + '_binary_' + str(feature_id),
            assetId=bin_asset_path,
            region=region,
            scale=scale,
            crs=crs,
            maxPixels=1e13,
            pyramidingPolicy={'.default': 'MODE'}
        ).start()
        print(f'  Submitted binary: {bin_asset_path}')

        if export_probability:
            prob_stack = pipeline.aggregate_probability_mc(binary_mc_ic, ds, year_list_feat)
            prob_stack = (prob_stack
                          .clip(region)
                          .set('feature_id', feature_id,
                               'years_processed', years_processed,
                               'n_years', len(years_processed))
                          .updateMask(prob_stack.mask().gt(0)))

            prob_asset_path = (output_base_path + probability_dir + '/'
                               + ds + '_probability_' + out_version + '/'
                               + str(feature_id))
            ee.batch.Export.image.toAsset(
                image=prob_stack,
                description=ds + '_probability_' + str(feature_id),
                assetId=prob_asset_path,
                region=region,
                scale=scale,
                crs=crs,
                maxPixels=1e13
            ).start()
            print(f'  Submitted probability: {prob_asset_path}')

    print(f'  All tasks submitted for feature {feature_id} '
          f'({len(years_processed)} years: {years_processed[0]}–{years_processed[-1]})')

print('Done. Check GEE task manager for results.')
