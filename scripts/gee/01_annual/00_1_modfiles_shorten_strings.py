"""

DESCRIPTION: Shorten predictor names in GCS decision-tree string files.

AUTHOR: Kathleen Orndahl
DATE: 2026-07-01

NOTES:

- Reads each tree file from 'final_trees/', applies a predictor name shortening map, and
  writes the shortened file to 'final_trees_shortened/' in the same GCS bucket.

- Three independent shortening steps are applied, together reducing file sizes by ~28%
  (empirically estimated on mc=1):

  1. Predictor name shortening: e.g. spectral_NDVI_annualMean -> s_NDVI_AM,
     ecoregion_CanadianMiddleArcticTundra -> er_CMAT  (~17% reduction).
     The same shorten() mapping is used in run_annual_biomass_tree_shorten.py to rename
     predictor image bands before classification. Both scripts must stay in sync.

  2. Threshold precision truncation: 15dp -> 4dp (e.g. 6413.72030450659 -> 6413.7203).
     Safe — spectral reflectance thresholds are integer-scaled so fractional digits beyond
     the integer are already meaningless; lat/lon at 4dp is still ~11 m precision.
     Produces identical classifications on all real input data. (~9% reduction)

  3. Trailing whitespace stripping: R's ranger output appends 2 spaces to ~50% of lines.
     (~1.2% reduction)

  NOTE: Each split-node line also contains two literal '0 0' columns (node sample count
  and deviance) that are always exactly zero in ranger random forest output. Removing them
  would save an additional ~5%, but requires testing with GEE's tree parser to confirm
  the format is accepted. Not implemented here.

- Idempotent: blobs that already exist in final_trees_shortened/ are skipped.

- Processes 4 model combinations x 100 MC = 400 files total.

- Run from scripts/gee/ with: python annual/shorten_tree_strings.py

PREREQUISITES:
  google-cloud-storage Python package installed.
  Application Default Credentials configured (gcloud auth application-default login).

"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import re
from google.cloud import storage
from utils import misc

# =========================
# 1. CONFIGURATION ========
# =========================

BUCKET_NAME = 'arctic_biomass_mapping_models'
MODEL_VERSION = 'v20240514'
MC_LIST = list(range(1, 101))
THRESHOLD_DP = 4

# =========================
# 2. SHORTENING MAP =======
# =========================

# Regex matching split-node lines:  {spaces}{node}) {predictor} {< or >} {threshold} ...
# Group 1: prefix (spaces + node id + ") ")
# Group 2: predictor name
# Group 3: comparison operator prefix (" <" or " >")
_SPLIT_PATTERN = re.compile(r'^(\s+\d+\)\s+)(\S+)(\s+[<>])', re.MULTILINE)

# Regex matching threshold values with more digits than THRESHOLD_DP.
# Matches any floating-point number (possibly negative) with >THRESHOLD_DP decimal digits,
# positioned immediately after the comparison operator and space on a split-node line.
_THRESHOLD_PATTERN = re.compile(
    r'^(\s+\d+\)\s+\S+\s+[<>]\s+)(-?\d+\.\d{' + str(THRESHOLD_DP + 1) + r',})',
    re.MULTILINE
)

def shorten_threshold_precision(content):
    """Truncate threshold values in tree split lines to THRESHOLD_DP decimal places."""
    def _truncate(m):
        return m.group(1) + f'{float(m.group(2)):.{THRESHOLD_DP}f}'
    return _THRESHOLD_PATTERN.sub(_truncate, content)


def strip_trailing_whitespace(content):
    """Remove trailing spaces/tabs from every line (R's ranger appends 2 spaces per line)."""
    return re.sub(r'[ \t]+$', '', content, flags=re.MULTILINE)


def shorten_tree_content(content):
    """Apply all three shortening steps to a tree file string."""
    def replace_name(m):
        return m.group(1) + misc.shortenPredName(m.group(2)) + m.group(3)
    content = _SPLIT_PATTERN.sub(replace_name, content)       # step 1: predictor names
    content = shorten_threshold_precision(content)             # step 2: threshold precision
    content = strip_trailing_whitespace(content)               # step 3: trailing whitespace
    return content


# =========================
# 3. MAIN =================
# =========================

def main():

    client = storage.Client()
    bucket = client.bucket(BUCKET_NAME)

    model_combinations = [
        ('binary',     'total'),
        ('binary',     'woody'),
        ('continuous', 'total'),
        ('continuous', 'woody'),
    ]

    total = 0
    skipped = 0
    processed = 0
    errors = 0

    for model_type, ds_type in model_combinations:
        print(f'\n--- {ds_type}_{model_type} ---')

        for mc in MC_LIST:

            fname = (f'{ds_type}_{model_type}_formatted_gee_short_{MODEL_VERSION}_{mc}.txt')
            src_path = f'{MODEL_VERSION}/mc/{model_type}/final_trees/{fname}'
            dst_path = f'{MODEL_VERSION}/mc/{model_type}/final_trees_shortened/{fname}'
            total += 1

            # Skip if destination already exists (safe to re-run)
            dst_blob = bucket.blob(dst_path)
            if dst_blob.exists():
                skipped += 1
                if mc == 1:
                    print(f'  mc={mc:03d}: skipping (exists)')
                continue

            try:
                src_blob = bucket.blob(src_path)
                content = src_blob.download_as_text(encoding='utf-8')
                shortened = shorten_tree_content(content)

                dst_blob.upload_from_string(shortened.encode('utf-8'), content_type='text/plain; charset=utf-8')

                orig_kb  = len(content.encode('utf-8')) / 1024
                short_kb = len(shortened.encode('utf-8')) / 1024
                saved_pct = (1 - short_kb / orig_kb) * 100
                print(f'  mc={mc:03d}: {orig_kb:.0f} KB -> {short_kb:.0f} KB  (-{saved_pct:.1f}%)')
                processed += 1

            except Exception as e:
                print(f'  mc={mc:03d}: ERROR — {e}')
                errors += 1

    print(f'\n=== Done ===')
    print(f'  Processed: {processed}')
    print(f'  Skipped (already exist): {skipped}')
    print(f'  Errors: {errors}')
    print(f'  Total files: {total}')


if __name__ == '__main__':
    main()
