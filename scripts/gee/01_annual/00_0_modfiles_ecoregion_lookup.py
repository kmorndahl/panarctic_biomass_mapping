"""

DESCRIPTION: Generate a lookup table mapping ecoregion predictor band names to short codes
             and save it to utils/params.py.

AUTHOR: Kathleen Orndahl
DATE: 2026-07-01

NOTES:

- Loads the ecoregions_img asset from GEE, sorts band names alphabetically, and assigns
  sequential codes er_001, er_002, ... to each ecoregion.

- Writes the resulting _ECOREGION_SHORT dict to utils/params.py so it can be imported
  by shorten_tree_strings.py and run_annual_biomass_tree_shorten.py. Re-running this
  script overwrites the existing dict in params.py.

- Run this script once whenever the ecoregions_img asset changes (band names added or
  removed). After running, re-run shorten_tree_strings.py to regenerate the shortened
  tree files.

- Run from scripts/gee/ with: python annual/ecoregion_name_lookup.py

PREREQUISITES:
  earthengine-api Python package installed.
  Authenticated via gcloud auth application-default login or ee.Authenticate().

"""

# =========================
# 0. PATH SETUP ===========
# =========================

import sys, os, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import ee
from utils import params

try:
    ee.Initialize(project=params.ee_project)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=params.ee_project)

# =========================
# 1. LOAD ECOREGION BANDS =
# =========================

band_names = sorted(
    ee.Image('projects/arctic-biomass-mapping/assets/predictors/ecoregions_img')
    .bandNames()
    .getInfo()
)

print(f'Found {len(band_names)} ecoregion bands.\n')

# =========================
# 2. BUILD LOOKUP TABLE ===
# =========================

mapping = {name: f'er_{i + 1:03d}' for i, name in enumerate(band_names)}

print(f'  {"Full name":<65s}  Short code')
print(f'  {"-" * 65}  ----------')
for name, code in mapping.items():
    print(f'  {name:<65s}  {code}')

# =========================
# 3. WRITE TO params.py ===
# =========================

lines = ['_ECOREGION_SHORT = {']
for name, code in mapping.items():
    lines.append(f"    '{name}': '{code}',")
lines.append('}')
new_block = '\n'.join(lines)

params_path = os.path.join(os.path.dirname(__file__), '..', 'utils', 'params.py')
with open(params_path, 'r') as f:
    content = f.read()

existing = re.search(r'^_ECOREGION_SHORT\s*=\s*\{[^}]*\}', content, re.MULTILINE)
if existing:
    content = content[:existing.start()] + new_block + content[existing.end():]
else:
    content = content.rstrip() + '\n\n' + new_block + '\n'

with open(params_path, 'w') as f:
    f.write(content)

print(f'\nSaved _ECOREGION_SHORT ({len(mapping)} entries) to utils/params.py.')
