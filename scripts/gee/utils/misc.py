import ee
import time
from re import sub
from utils import params

# FUNCTION: calcAnnualSpectral
# USE: Calculate annual summaries (mean, median, min, max, range) for modeled reflectance data
# PARAMETERS:
#  spectral = modeled seasonal reflectance data
#  band = spectral band/index
# NOTES: Change rates (formerly computed here per-band) are now computed separately
#        for all bands at once via calcAllChangeRates() — see that function for why.
# AUTHOR: Kathleen Orndahl

def calcAnnualSpectral(spectral, band):

    regex = '.*' + band + '.*';
    annual_spectral = spectral.select(regex);

    # Get annual summaries — mean, median, min, max in a single reduce() call instead of three.
    # combine() output band order is fixed based on input order: [mean, median, min, max], so we rename positionally
    # rather than relying on inferred band names.
    combined_reducer = (ee.Reducer.mean()
                         .combine(ee.Reducer.median(), sharedInputs=True)
                         .combine(ee.Reducer.minMax(), sharedInputs=True))
    summaries = annual_spectral.reduce(combined_reducer).rename([
        'spectral_' + band + '_annualMean',
        'spectral_' + band + '_annualMedian',
        'spectral_' + band + '_min',
        'spectral_' + band + '_max'
    ])
    mean = summaries.select('spectral_' + band + '_annualMean')
    median = summaries.select('spectral_' + band + '_annualMedian')
    annual_range = summaries.select('spectral_' + band + '_max') \
                            .subtract(summaries.select('spectral_' + band + '_min')) \
                            .rename('spectral_' + band + '_annualRange');

    # Change rates are computed across all bands at once in calcAllChangeRates() —
    # not per-band here — to keep the GEE graph compact.
    return annual_spectral.addBands(mean) \
                          .addBands(median) \
                          .addBands(annual_range)

# FUNCTION: calcChangeRate
# USE: Calculate rate of modeled reflectance change between seasons
# PARAMETERS:
#  annual_spectral = modeled seasonal reflectance data
#  band = spectral band/index
#  doys = seasonal DOY data
#  start = start season
#  end = end season
# NOTES: Called within calcAnnualSpectral, parameters populated there
# AUTHOR: Katie Orndahl

def calcChangeRate(annual_spectral, band, doys, start, end):

    doy_start = doys.select('doy_' + start)
    doy_end = doys.select('doy_' + end)

    doy_diff = doy_end.subtract(doy_start)

    spectral_start = annual_spectral.select('spectral_' + band + '_' + start)
    spectral_end = annual_spectral.select('spectral_' + band + '_' + end)

    spectral_diff = spectral_end.subtract(spectral_start)

    start = sub("[a-z]", '', start[0].upper() + start[1:])
    end = sub("[a-z]", '',  end[0].upper() + end[1:])

    return spectral_diff.divide(doy_diff).rename('spectral_' + band + '_change' + start + end)

# FUNCTION: calcAllChangeRates
# USE: Compute seasonal change rates for ALL spectral bands in four batched operations.
#      Replaces the 44 per-band calls that calcChangeRate() previously made inside
#      calcAnnualSpectral() — each .subtract()/.divide() here operates on an 11-band image
#      (one node in the GEE graph) rather than producing 11 separate single-band nodes.
# PARAMETERS:
#  spectral = modeled seasonal reflectance (all bands, all seasons) — same image passed to
#             calcAnnualSpectral(); must contain bands named spectral_<band>_<season>
#  doys     = seasonal DOY mosaic (ccdc_data['seasonal_doy'])
# AUTHOR: Kathleen Orndahl

def calcAllChangeRates(spectral, doys):

    # Season pairs with their abbreviated suffix, matching the naming convention used in
    # calcChangeRate() (sub("[a-z]", '') on the capitalised season name).
    season_pairs = [
        ('startSnowfree', 'earlySummer', 'SSES'),
        ('earlySummer',   'peakSummer',  'ESPS'),
        ('peakSummer',    'lateSummer',  'PSLS'),
        ('lateSummer',    'endSnowfree', 'LSES'),
    ]

    change_rate_images = []
    for start, end, suffix in season_pairs:

        # DOY difference between the two seasons (1-band image; broadcast across all spectral bands)
        doy_diff = doys.select('doy_' + end).subtract(doys.select('doy_' + start))

        # Select all spectral bands for each season at once — ONE subtract node for all 11 bands
        start_bands = spectral.select('.*_' + start)
        end_bands   = spectral.select('.*_' + end)
        change_rate = end_bands.subtract(start_bands).divide(doy_diff)

        # Rename: e.g. 'spectral_blue_earlySummer' → 'spectral_blue_changeSSES'
        change_rate = change_rate.regexpRename('_' + end + '$', '_change' + suffix)
        change_rate_images.append(change_rate)

    # Stack the four season-pair results (4 × 11 = 44 bands total)
    # Looping like this avoids annoying ee.ImageCollection().toBands() behavior
    # that prefixes image IDs
    result = change_rate_images[0]
    for img in change_rate_images[1:]:
        result = result.addBands(img)
    return result

# FUNCTION: createCollectionIfNotExists
# USE: Checks if a GEE ImageCollection exists. If not, it creates an empty one.
# PARAMETERS:
#  asset_id = path to check for existing ImageCollection
# AUTHOR: Katie Orndahl
def createCollectionIfNotExists(asset_id):

    try:
        # Attempt to get the asset metadata
        ee.data.getAsset(asset_id)
        print(f"Asset already exists: {asset_id}")

    except ee.EEException as e:
        # If it throws an EEException, it likely doesn't exist
        if "not found" in str(e).lower() or "not exist" in str(e).lower():
            print(f"Asset not found. Creating empty ImageCollection at: {asset_id}")

            # Create the empty ImageCollection asset
            ee.data.createAsset(
                {'type': 'IMAGE_COLLECTION'},
                asset_id
            )
        else:
            # If it's a different Earth Engine error, raise it
            raise e

# FUNCTION: createDirectoryIfNotExists
# USE: Checks if a GEE directory exists. If not, it creates an empty one.
# PARAMETERS:
#  asset_id = path to check for existing directory
# AUTHOR: Katie Orndahl
def createDirectoryIfNotExists(asset_id):

    try:
        # Attempt to get the asset metadata
        ee.data.getAsset(asset_id)
        print(f"Asset already exists: {asset_id}")

    except ee.EEException as e:
        # If it throws an EEException, it likely doesn't exist
        if "not found" in str(e).lower() or "not exist" in str(e).lower():
            print(f"Asset not found. Creating empty directory at: {asset_id}")

            # Create the empty ImageCollection asset
            ee.data.createAsset(
                {'type': 'FOLDER'},
                asset_id
            )
        else:
            # If it's a different Earth Engine error, raise it
            raise e
        
def shortenPredName(name):
    """
    Apply predictor name shortening to a single name string.

    Substitution order matters: prefix replacements come before suffix replacements so that
    e.g. 'spectral_NDVI_annualMean' -> 's_NDVI_annualMean' -> 's_NDVI_AM'.
    Ecoregion bands are handled first (early return) to avoid interference with other rules.
    """
    n = name

    # ----- Ecoregion bands (early return — sequential lookup from _ECOREGION_SHORT) -----
    # Lookup table generated by annual/ecoregion_name_lookup.py and stored in utils/params.py.
    # Run that script before this one if the ecoregions_img asset band names have changed.
    if n.startswith('ecoregion_'):
        return params._ECOREGION_SHORT[n]

    # ----- WTE class bands (full prefix before generic spectral_ replacement) -----
    n = sub(r'world_terrestrial_ecosystems_X', 'wte_X', n)

    # ----- Spectral / texture prefixes -----
    n = n.replace('spectral_', 's_')
    n = n.replace('texture_NDVI_peakSummer', 't_NDVI_PS')  # full name match

    # ----- Latitudinal zone bands -----
    n = n.replace('zone_3_LowArctic',  'z_LA')
    n = n.replace('zone_3_OroArctic',  'z_OA')
    n = n.replace('zone_3_SubArctic',  'z_SA')
    n = n.replace('zone_3_HighArctic', 'z_HA')

    # ----- Topographic predictors -----
    n = n.replace('topo_eastness',  'tp_east')
    n = n.replace('topo_northness', 'tp_north')
    n = n.replace('topo_slope',     'tp_slp')
    n = n.replace('topo_tpi',       'tp_tpi')
    n = n.replace('topo_dem',       'tp_dem')
    n = n.replace('topo_cti',       'tp_cti')

    # ----- Other predictors -----
    n = n.replace('permafrost_index', 'pf_idx')
    n = n.replace('trees_presence',   't_pres')
    n = n.replace('latitude',         'lat')
    n = n.replace('longitude',        'lon')

    # ----- Season / stat suffixes (applied after prefix substitutions) -----
    n = n.replace('_annualMedian',  '_AMD')
    n = n.replace('_annualMean',    '_AM')
    n = n.replace('_annualRange',   '_AR')
    n = n.replace('_startSnowfree', '_SS')
    n = n.replace('_earlySummer',   '_ES')
    n = n.replace('_peakSummer',    '_PS')
    n = n.replace('_lateSummer',    '_LS')
    n = n.replace('_endSnowfree',   '_EN')
    n = n.replace('_changePSLS',    '_cPSLS')
    n = n.replace('_changeSSES',    '_cSSES')
    n = n.replace('_changeESPS',    '_cESPS')
    n = n.replace('_changeLSES',    '_cLSES')

    return n


# FUNCTION: wait_for_task_slots
# USE: Block task submission until the number of active GEE export tasks drops below a
#      threshold, preventing the per-project ~5000 concurrent task limit from being exceeded.
#      A time-gate skips the listOperations API call when checked recently, so overhead is
#      at most one call per min_check_interval seconds regardless of submission rate.
# PARAMETERS:
#  project            = GEE project ID string (e.g. params.ee_project_production)
#  max_active         = pause when active task count reaches this threshold (default 4500)
#  poll_interval      = seconds to wait between rechecks when the queue is full
#  min_check_interval = minimum seconds between listOperations calls in the happy path
# AUTHOR: Kathleen Orndahl

_last_task_check = 0.0  # module-level timestamp; persists across calls within a session

def wait_for_task_slots(project, max_active=4500, poll_interval=60, min_check_interval=30):
    global _last_task_check
    project_path = f'projects/{project}'
    if time.time() - _last_task_check < min_check_interval:
        return
    while True:
        ops = ee.data.listOperations(project_path)
        n_active = sum(1 for op in ops if not op.get('done', False))
        _last_task_check = time.time()
        if n_active < max_active:
            break
        print(f'  Task queue: {n_active} active (limit {max_active}). Waiting {poll_interval}s...')
        time.sleep(poll_interval)