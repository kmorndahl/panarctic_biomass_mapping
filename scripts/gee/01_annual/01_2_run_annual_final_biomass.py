"""

DESCRIPTION: Apply binary presence/absence masking and ecological constraint to produce final biomass maps.

AUTHOR: Kathleen Orndahl
DATE: 2026-07-06

NOTES:

- Run AFTER all 01_1_run_annual_binary_continuous.py exports have completed in GEE. All four
  input collections (total_continuous, total_binary, woody_continuous, woody_binary) must exist
  before this script runs.

WHY THIS SCRIPT EXISTS:

  01_1_run_annual_binary_continuous.py exports continuous and binary models as SEPARATE tasks to
  stay within GEE's server-side encoded-object size limit per tile. This script combines them:

    1. total_biomass  = total_continuous  × total_presence_absence
    2. woody_biomass  = woody_continuous  × total_presence_absence × woody_presence_absence
    3. Ecological constraint: total_biomass_median = max(total_biomass_median, woody_biomass_median)
       Applied to the median (p50) band only. Uncertainty bounds (lwr, upr) reflect the
       raw MC distribution and are not constrained.

  Since all inputs are already-exported asset leaf nodes, GEE evaluates this as a trivial
  expression graph — no decision-tree strings are loaded and exports complete in seconds per tile.

LOOP STRUCTURE:

  Driven by the completed total_continuous ImageCollection (one image per feature, all years
  as bands). For each feature, a Python year-loop selects per-year bands, applies masking and
  the ecological constraint, and stacks all years into a single output image. Skips any feature
  where one of the other three inputs is not yet available — safe to re-run as more Phase 2
  tasks complete.

OUTPUT (one image per feature, all years as bands):

  final_biomass/total_biomass_{version}/ — 3 bands per year:
    total_biomass_{year}, total_biomass_lwr_{year}, total_biomass_upr_{year}

  final_biomass/woody_biomass_{version}/ — 3 bands per year:
    woody_biomass_{year}, woody_biomass_lwr_{year}, woody_biomass_upr_{year}

"""

# =========================
# 0. PATH SETUP ===========
# =========================

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import ee
from utils import params

try:
    ee.Initialize(project=params.ee_project)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=params.ee_project)

# =========================
# 1. SET UP ===============
# =========================

# 1.0 ----- PARAMETERS -----

# Reprojection options
crs = params.crs
scale = params.scale

# Versions — must match out_version from 01_1_run_annual_binary_continuous.py
out_version = '20260705'
output_dir = 'annual_comparison'
output_base_path = 'projects/arctic-biomass-change/assets/' + output_dir + '/'

# Input directories (from 01_1_run_annual_binary_continuous.py)
continuous_dir = 'continuous_TEST'
binary_dir = 'binary'

# Output directory (produced by this script)
final_biomass_dir = 'final_biomass'

# 1.1 ----- MODULES -----

from utils import misc as fun_misc

# 1.2 ----- CREATE OUTPUT COLLECTIONS -----

fun_misc.createDirectoryIfNotExists(output_base_path + final_biomass_dir)
fun_misc.createCollectionIfNotExists(
    output_base_path + final_biomass_dir + '/total_biomass_' + out_version
)
fun_misc.createCollectionIfNotExists(
    output_base_path + final_biomass_dir + '/woody_biomass_' + out_version
)

# 1.3 ----- LOAD INPUT IMAGE COLLECTIONS -----

# Load all four input collections produced by 01_1_run_annual_binary_continuous.py.
# Use aggregate_array() rather than toList().getInfo() — images have up to 120 bands
# each and serializing full image objects for 1,000 features can exceed GEE user memory.
total_cont_ic = ee.ImageCollection(
    output_base_path + continuous_dir + '/total_continuous_' + out_version
)
woody_cont_ic = ee.ImageCollection(
    output_base_path + continuous_dir + '/woody_continuous_' + out_version
)
total_bin_ic = ee.ImageCollection(
    output_base_path + binary_dir + '/total_binary_' + out_version
)
woody_bin_ic = ee.ImageCollection(
    output_base_path + binary_dir + '/woody_binary_' + out_version
)

# Fetch only the scalar/list properties needed — one round-trip per collection
total_cont_fids  = total_cont_ic.aggregate_array('feature_id').getInfo()
total_cont_ids   = total_cont_ic.aggregate_array('system:id').getInfo()
total_cont_years = total_cont_ic.aggregate_array('years_processed').getInfo()

woody_cont_fids  = woody_cont_ic.aggregate_array('feature_id').getInfo()
woody_cont_ids   = woody_cont_ic.aggregate_array('system:id').getInfo()

total_bin_fids   = total_bin_ic.aggregate_array('feature_id').getInfo()
total_bin_ids    = total_bin_ic.aggregate_array('system:id').getInfo()

woody_bin_fids   = woody_bin_ic.aggregate_array('feature_id').getInfo()
woody_bin_ids    = woody_bin_ic.aggregate_array('system:id').getInfo()

# Build feature_id → asset id lookups for the three non-driver collections
woody_cont_by_fid = dict(zip(woody_cont_fids, woody_cont_ids))
total_bin_by_fid  = dict(zip(total_bin_fids,  total_bin_ids))
woody_bin_by_fid  = dict(zip(woody_bin_fids,  woody_bin_ids))

# year list per feature (years_processed was stored as a list property in Phase 2)
years_by_fid = {fid: [int(y) for y in yrs]
                for fid, yrs in zip(total_cont_fids, total_cont_years)}

total_cont_by_fid = dict(zip(total_cont_fids, total_cont_ids))

# 1.4 ----- REPORT -----

n_complete = sum(
    1 for fid in total_cont_fids
    if fid in woody_cont_by_fid and fid in total_bin_by_fid and fid in woody_bin_by_fid
)

print(f'Total continuous features available: {len(total_cont_fids)}')
print(f'Features with all 4 inputs complete: {n_complete}')

# ==============================
# 2. MAIN LOOP — FEATURES ======
# ==============================

for feature_id in total_cont_fids:

    # 2.1 ----- CHECK COMPLETENESS -----

    if feature_id not in woody_cont_by_fid:
        print(f'Skipping {feature_id}: woody_continuous not yet exported')
        continue
    if feature_id not in total_bin_by_fid:
        print(f'Skipping {feature_id}: total_binary not yet exported')
        continue
    if feature_id not in woody_bin_by_fid:
        print(f'Skipping {feature_id}: woody_binary not yet exported')
        continue

    year_list = years_by_fid[feature_id]
    print(f'Processing feature {feature_id} ({len(year_list)} years: {year_list[0]}–{year_list[-1]})')

    # 2.2 ----- LOAD INPUTS AS ASSET LEAF NODES -----

    # All four are ee.Image asset references — trivial leaf nodes; no tree strings are loaded.
    total_cont = ee.Image(total_cont_by_fid[feature_id])
    woody_cont = ee.Image(woody_cont_by_fid[feature_id])
    total_pa   = ee.Image(total_bin_by_fid[feature_id])
    woody_pa   = ee.Image(woody_bin_by_fid[feature_id])
    region = total_cont.geometry()

    # 2.3 ----- APPLY MASKING AND ECOLOGICAL CONSTRAINT PER YEAR -----

    # Python year loop selects per-year bands, applies masking, and assembles the
    # constrained output. All select/multiply/max calls are lazy (no data fetched).
    total_year_bands = []
    woody_year_bands = []

    for year in year_list:

        y = str(year)

        # Per-year band selection
        t_cont_y = total_cont.select([
            'total_biomass_'     + y,
            'total_biomass_lwr_' + y,
            'total_biomass_upr_' + y
        ])
        w_cont_y = woody_cont.select([
            'woody_biomass_'     + y,
            'woody_biomass_lwr_' + y,
            'woody_biomass_upr_' + y
        ])
        t_pa_y = total_pa.select('total_presence_absence_' + y)   # 1-band
        w_pa_y = woody_pa.select('woody_presence_absence_' + y)   # 1-band

        # Masking — multiplying by a 1-band presence/absence image broadcasts
        # across all 3 bands of the continuous image
        t_masked_y = t_cont_y.multiply(t_pa_y)
        w_masked_y = w_cont_y.multiply(t_pa_y).multiply(w_pa_y)

        # Ecological constraint: enforce total >= woody at the median.
        # Uncertainty bounds are left unchanged — they reflect the raw MC distribution.
        t_median_constrained = (t_masked_y.select('total_biomass_' + y)
                                .max(w_masked_y.select('woody_biomass_' + y)))

        t_final_y = (t_median_constrained
                     .addBands(t_masked_y.select('total_biomass_lwr_' + y))
                     .addBands(t_masked_y.select('total_biomass_upr_' + y)))

        total_year_bands.append(t_final_y)
        woody_year_bands.append(w_masked_y)

    # 2.4 ----- STACK ALL YEARS -----

    total_final = (ee.ImageCollection(total_year_bands)
                   .toBands()
                   .regexpRename('^[0-9]+_', '')
                   .round().uint16()
                   .set('feature_id', feature_id,
                        'years_processed', year_list,
                        'n_years', len(year_list)))

    woody_final = (ee.ImageCollection(woody_year_bands)
                   .toBands()
                   .regexpRename('^[0-9]+_', '')
                   .round().uint16()
                   .set('feature_id', feature_id,
                        'years_processed', year_list,
                        'n_years', len(year_list)))

    # 2.5 ----- EXPORT FINAL BIOMASS MAPS -----

    for ds, img in [('total', total_final), ('woody', woody_final)]:

        asset_path = (output_base_path + final_biomass_dir + '/' + ds + '_biomass_'
                      + out_version + '/' + str(feature_id))

        ee.batch.Export.image.toAsset(
            image=img,
            description=ds + '_biomass_' + str(feature_id),
            assetId=asset_path,
            region=region,
            scale=scale,
            crs=crs,
            maxPixels=1e13
        ).start()

        print(f'  Submitted {ds}_biomass: {asset_path}')

print('All final biomass tasks submitted.')
