"""

DESCRIPTION: Run reflectance calculation pipeline — aggregated-years variant

AUTHOR: Kathleen Orndahl
DATE: 2026-07-05

NOTES:

- Restructured from 01_0_run_annual_reflectance.py to export ONE image per feature containing
  all years as bands, instead of one image per feature-year. Reduces task count from
  N_features × N_years (40,000) to N_features (1,000), allowing 40× more work per concurrent
  GEE task slot.

- Phase 2 (01_1_run_annual_binary_continuous.py) reads these multi-year assets and
  extracts per-year bands using band name suffix matching (e.g. '.*_1984').

- Band naming convention: <original_band_name>_<year>
    e.g. spectral_blue_startSnowfree_1984, spectral_NDVI_peakSummer_ND_1984, ...

LOOP STRUCTURE:

  feature (outer) × year (inner, Python)

  Year-specific lazy images (refl, refl_tc) are built once per year-within-feature and
  renamed with a year suffix before stacking. GEE evaluates the full computation only once
  per tile during the single export task, not once per year-task.

TASK COUNT: N_features = 1,000 (vs. N_features × N_years = 40,000 in the per-year script)

"""

# =========================
# 0. PATH SETUP ===========
# =========================

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import ee
from utils import params
import time

try:
    ee.Initialize(project=params.ee_project_production)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=params.ee_project_production)

# =========================
# 1. SET UP ===============
# =========================

# 1.0 ----- PARAMETERS -----

# Feature collection
features_fc_path = 'projects/arctic-biomass-change/assets/ROIs/landsat_sample_points_n50000'
subset = True   # Whether to randomly subset feature collection
n_feat = 1000   # If subset=True, number of features to select
feat_buffer = 120  # Buffer around features in meters; set to 0 for no buffer

# Property name used to build output asset IDs
feature_id_property = 'sample_id'

# Year list — all years processed together per feature in a single export task
year_list = list(range(1984, 2024))

# CCDC options
ccdc_version = 'v20240207'
ccdc_start_year = 1984
ccdc_end_year = 2023

# CCDC processing options
extrapolate_max_days = 120
segment_find_strategy = 'next'

# Topographic correction parameters
topo_seasons = ['startSnowfree', 'earlySummer', 'peakSummer', 'lateSummer', 'endSnowfree']
bands_topo_corr = ['spectral_blue_', 'spectral_green_', 'spectral_red_', 'spectral_NIR_',
                   'spectral_SWIR1_', 'spectral_SWIR2_', 'spectral_EVI2b_']
topo_buffer = 9000  # Meters; feature is buffered by this amount for the topo correction fit

# Reprojection options
crs = params.crs
scale = params.scale

# Output options
out_version = '20260705'
output_dir = 'annual_comparison'
output_base_path = 'projects/arctic-biomass-change/assets/' + output_dir + '/'
product_category = 'modeled_reflectance'

# Processing options
tile_scale = 1

# 1.1 ----- MODULES -----

from utils import pipeline
from utils import misc as fun_misc

# 1.2 ----- YEAR-INDEPENDENT DATA -----

print('Loading static predictors and CCDC data...')
static = pipeline.load_static_predictors()
ccdc_data = pipeline.load_ccdc_data(ccdc_version, ccdc_start_year, ccdc_end_year)

# 1.3 ----- CREATE OUTPUT COLLECTION -----

fun_misc.createDirectoryIfNotExists(output_base_path + product_category)
fun_misc.createCollectionIfNotExists(
    output_base_path + product_category + '/' + product_category + '_' + out_version
)

# 1.4 ----- LOAD FEATURE COLLECTION -----

features_fc = ee.FeatureCollection(features_fc_path)

if subset:
    features_fc = (features_fc
        .randomColumn(columnName='random', seed=42)
        .sort('random')
        .limit(n_feat))

# Fetch all features client-side in one getInfo() call
feature_id_list = features_fc.aggregate_array(feature_id_property).getInfo()
features_list = features_fc.toList(len(feature_id_list)).getInfo()
features_by_id = {f['properties'][feature_id_property]: f for f in features_list}

# 1.5 ----- REPORT -----

print(f'Total features to process: {len(feature_id_list)}')
print(f'Total years per feature: {len(year_list)}  ({year_list[0]}–{year_list[-1]})')
print(f'Total tasks to submit: {len(feature_id_list)}')

feature_id_list = ['S_12107', 'S_489', 'S_26776', 'S_46791', 'S_28375', 'S_37513', 'S_4690', 'S_44283', 'S_15732', 'S_15215', 'S_10123', 'S_13432', 'S_42925', 'S_2587', 'S_39840', 'S_3570', 'S_14015', 'S_35429', 'S_41214', 'S_27523', 'S_36499', 'S_15625', 'S_20561', 'S_42267', 'S_45035', 'S_41867', 'S_23171', 'S_48668', 'S_30029', 'S_19883', 'S_330', 'S_30846', 'S_33297', 'S_10473', 'S_5709', 'S_7671', 'S_27026', 'S_48872', 'S_42940', 'S_3169', 'S_21321', 'S_40968', 'S_25293', 'S_20432', 'S_19603', 'S_14054', 'S_2059', 'S_33754', 'S_42032', 'S_16054', 'S_46133', 'S_32037', 'S_7616', 'S_37854', 'S_17192', 'S_14024', 'S_12450', 'S_26924', 'S_4140', 'S_41213', 'S_40388', 'S_23342', 'S_2762', 'S_26481', 'S_32050', 'S_16009', 'S_24758', 'S_20163', 'S_25869', 'S_31380', 'S_45842', 'S_40513', 'S_21782', 'S_34277', 'S_27375', 'S_20940', 'S_21038', 'S_10716', 'S_894', 'S_38519', 'S_5514', 'S_1933', 'S_44795', 'S_16452', 'S_6883', 'S_10981']

# ==============================
# 2. OUTER LOOP — FEATURES =====
# ==============================

for feature_id in feature_id_list:

    start_time = time.perf_counter()

    print(f'===== FEATURE: {feature_id} =====')

    # -------------------------------------------------------------------------
    # 2.1 FEATURE GEOMETRY
    # -------------------------------------------------------------------------

    feature = ee.Feature(features_by_id[feature_id])
    if feat_buffer > 0:
        feature = feature.buffer(feat_buffer)

    roi_topo = ee.FeatureCollection(feature.buffer(topo_buffer))

    # -------------------------------------------------------------------------
    # 2.2 SEASONAL REFLECTANCE — ALL YEARS (lazy, no computation triggered yet)
    # -------------------------------------------------------------------------

    refl_by_year = {}
    for year in year_list:
        refl_by_year[year] = pipeline.compute_seasonal_reflectance(
            year, ccdc_data, extrapolate_max_days, segment_find_strategy
        )

    # -------------------------------------------------------------------------
    # 2.3 TOPOGRAPHIC CORRECTION COEFFICIENTS
    # eecu: one batched getInfo() per feature; c values embedded as float constants
    #       → no lazy linearFit in the export graph → 1–10 EECU/task.
    # time: deferred to export graph; fast submission but ~40 × 35 lazy ops per tile
    #       → causes per-tile memory/time overload → GEE task retries.
    # -------------------------------------------------------------------------

    c_values = pipeline.precompute_topo_c(
        refl_by_year, static['slp'], ccdc_data['ic'],
        topo_seasons, bands_topo_corr, roi_topo, crs, scale, tile_scale
    )

    # -------------------------------------------------------------------------
    # 2.4 INNER LOOP — YEARS
    # Apply topo correction and rename bands with year suffix before stacking.
    # -------------------------------------------------------------------------

    year_images = []

    for year in year_list:

        refl = refl_by_year[year]

        refl_tc = pipeline.apply_topo_correction_from_c(
            refl, static['slp'], ccdc_data['ic'],
            topo_seasons, bands_topo_corr, year, c_values
        )

        # Assemble: topo-corrected reflectance + normalized indices (not topo-corrected)
        refl_nd = refl.select('.*_ND.*|.*NBR.*')
        modeled_reflectance = refl_tc.addBands(refl_nd)

        # Append year suffix to all band names:
        #   spectral_blue_startSnowfree  →  spectral_blue_startSnowfree_1984
        year_suffix = '_' + str(year)
        year_refl = modeled_reflectance.rename(
            modeled_reflectance.bandNames().map(lambda b: ee.String(b).cat(year_suffix))
        )

        year_images.append(year_refl)

    # -------------------------------------------------------------------------
    # 2.5 STACK YEARS AND EXPORT
    # toBands() adds a numeric IC index prefix (e.g. '0_', '1_'); strip it so
    # final band names are clean: spectral_blue_startSnowfree_1984, etc.
    # -------------------------------------------------------------------------

    all_years_img = (ee.ImageCollection.fromImages(year_images)
                     .toBands()
                     .regexpRename('^[0-9]+_', ''))

    exp_img = (all_years_img
               .clip(feature.geometry())
               .set('feature_id', feature_id,
                    'year_start', year_list[0],
                    'year_end', year_list[-1],
                    'n_years', len(year_list))
               .updateMask(all_years_img.mask().gt(0)))

    asset_path = (output_base_path + product_category + '/'
                  + product_category + '_' + out_version + '/'
                  + str(feature_id))

    task = ee.batch.Export.image.toAsset(
        image=exp_img,
        description=product_category + '_' + str(feature_id),
        assetId=asset_path,
        region=feature.geometry(),
        scale=scale,
        crs=crs,
        maxPixels=1e13
    )
    task.start()

    print(f'  Submitted: {asset_path}')

    end_time = time.perf_counter()
    duration = end_time - start_time
    
    print(f"-> Iteration took {duration:.4f} seconds\n")

print('All tasks submitted.')
