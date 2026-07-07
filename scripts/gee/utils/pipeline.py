"""

DESCRIPTION: Shared pipeline functions for annual_sample and annual feature collection workflows

AUTHOR: Kathleen Orndahl
DATE: 2026-06-26

NOTES:

- All functions are pure: they take inputs and return ee objects, with no hardcoded asset paths
  except for fixed external assets (GCS model blobs, predictor catalogs) that are unlikely to change
- Functions return lazy ee.Image / ee.ImageCollection objects; no computation is triggered until
  an export or .getInfo() call is made by the calling script
- Import from scripts that have already inserted scripts/gee/ into sys.path:
    from utils import pipeline
- The Segments class and ccdc_tools must be accessible via the calling script's sys.path

"""

import ee

from utils import misc as fun_misc
from utils import refl as fun_refl
from utils import temporal_segmentation as temporal_seg
from ccdc_tools import dates as utils_dates
import math

Segments = temporal_seg.Segments


# =========================
# 1. DATA LOADERS =========
# =========================

# 1.0 ----- STATIC PREDICTORS -----

def load_static_predictors():
    """
    Load year-independent predictor data.

    Returns a dict of lazy ee.Image objects:
      'topographic' — multi-band MERIT DEM-based topographic stack (DEM, CTI, slope, TPI, eastness, northness);
                      gap-filled with a 3-pixel focal mean then blended with original
      'permafrost'  — Gruber (2012) permafrost zonation index
      'wte'         — World Terrestrial Ecosystems (encoded, gap-filled)
      'ecoregions'  — encoded ecoregion image (used in continuous model only)
      'zones'       — encoded Arctic latitudinal zone image
      'slp'         — Copernicus 30m slope (used for topographic correction in apply_topo_correction)
    """
    # Topographic predictors (MERIT DEM, 90m; gap-filled by focal mean then blended with original)
    dem = ee.Image("MERIT/DEM/v1_0_3").rename('topo_dem')
    cti = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/cti").mosaic().rename('topo_cti')
    slope_topo = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/slope").mosaic().rename('topo_slope')
    tpi = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/tpi").mosaic().rename('topo_tpi')
    eastness = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/eastness").mosaic().rename('topo_eastness')
    northness = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/northness").mosaic().rename('topo_northness')

    topographic = (dem.addBands(cti)
                      .addBands(slope_topo)
                      .addBands(tpi)
                      .addBands(eastness)
                      .addBands(northness))
    topographic_fill_gaps = topographic.focalMean(
        radius=dem.projection().nominalScale().multiply(3).divide(2),
        kernelType='square',
        units='meters',
        iterations=1
    )
    topographic = topographic_fill_gaps.blend(topographic)

    permafrost = ee.Image(
        'projects/arctic-biomass-mapping/assets/predictors/global_permafrost_zonation_index_gruber2012_filled_20km_v3'
    ).rename('permafrost_index')

    wte = ee.Image('projects/arctic-biomass-mapping/assets/predictors/wte_2020_filled_encoded')
    wte = wte.updateMask(wte.mask().gt(0))

    ecoregions = ee.Image('projects/arctic-biomass-mapping/assets/predictors/ecoregions_img')
    zones = ee.Image('projects/arctic-biomass-mapping/assets/predictors/zones_img')

    # Copernicus 30m slope — used as input to apply_topo_correction
    slp = ee.ImageCollection(
        'projects/arctic-biomass-mapping/assets/predictors/copernicus_30m_slope_gcs'
    ).mosaic()

    return {
        'topographic': topographic,
        'permafrost': permafrost,
        'wte': wte,
        'ecoregions': ecoregions,
        'zones': zones,
        'slp': slp
    }


# 1.1 ----- CCDC DATA -----

def load_ccdc_data(ccdc_version, start_year, end_year):
    """
    Load CCDC-derived data and topo correction illumination condition.

    Returns a dict:
      'ccdc_fit'    — mosaiced CCDC fit ImageCollection
      'seasonal_doy'— mosaiced seasonal DOY ImageCollection
      'doy_median'  — ee.Dictionary of median DOY per season (year-independent reference dates)
      'segments'    — Segments object for find_by_date (year-independent; date varies per year)
      'ic'          — illumination condition mosaic (year-independent; used in apply_topo_correction)
    """
    seasonal_doy = ee.ImageCollection(
        'projects/arctic-biomass-mapping/assets/seasonal_percentile_doys/seasonal_doys_tiles_'
        + str(start_year) + '_' + str(end_year) + '_' + ccdc_version
    ).mosaic()
    ccdc_fit = ee.ImageCollection(
        'projects/arctic-biomass-mapping/assets/CCDC_tiles/CCDC_C2_SR_tiles_'
        + str(start_year) + '_' + str(end_year) + '_' + ccdc_version
    ).mosaic()
    doy_median_fc = ee.FeatureCollection(
        'projects/arctic-biomass-mapping/assets/seasonal_percentile_doys/seasonal_doys_median_'
        + str(start_year) + '_' + str(end_year) + '_' + ccdc_version
    )
    ic = ee.ImageCollection(
        'projects/arctic-biomass-mapping/assets/topographic_correction/illumination_condition_'
        + str(start_year) + '_' + str(end_year) + '_' + ccdc_version
    ).mosaic()

    doy_median = doy_median_fc.first().toDictionary([
        'doy_earlySummer', 'doy_endSnowfree', 'doy_lateSummer', 'doy_peakSummer', 'doy_startSnowfree'
    ])
    segments = Segments(ccdc_fit, date_format=1)

    return {
        'ccdc_fit': ccdc_fit,
        'seasonal_doy': seasonal_doy,
        'doy_median': doy_median,
        'segments': segments,
        'ic': ic
    }


# ==================================
# 2. REFLECTANCE CALCULATION =======
# ==================================

# 2.0 ----- SEASONAL REFLECTANCE -----

def compute_seasonal_reflectance(year, ccdc_data, extrapolate_max_days, segment_find_strategy):
    """
    Compute CCDC-modeled seasonal reflectance for a given year.

    Port of 04_2 logic. Runs fitT() for each of the 5 seasons using the median DOY of each
    season as the reference date for segment selection.

    Parameters
    ----------
    year                 : int — calendar year to model reflectance for
    ccdc_data            : dict returned by load_ccdc_data()
    extrapolate_max_days : int — days to extrapolate beyond CCDC segment endpoints
    segment_find_strategy: str — 'next' for image exports; 'closest' for cal/val extraction

    Returns
    -------
    ee.Image — lazy multiband int16 image with bands spectral_<band>_<season>
    """
    seasonal_doy = ccdc_data['seasonal_doy']
    doy_median = ccdc_data['doy_median']
    segments = ccdc_data['segments']

    # Redefined per year because it captures 'year' by closure
    def convert_dates(key, value):
        frac_date = ee.Number(year).add(ee.Number(value).divide(365.25))
        return utils_dates.convertDate({'inputFormat': 1, 'inputDate': frac_date, 'outputFormat': 4})

    doy_median_frac = doy_median.map(convert_dates)

    start_snowfree_fit = fun_refl.fitT(
        seasonal_doy.select('doy_startSnowfree'), year,
        segments.find_by_date(ee.Date(doy_median_frac.get('doy_startSnowfree')), strategy=segment_find_strategy),
        extrapolate_max_days
    ).regexpRename('$', '_startSnowfree').regexpRename('^', 'spectral_')

    early_summer_fit = fun_refl.fitT(
        seasonal_doy.select('doy_earlySummer'), year,
        segments.find_by_date(ee.Date(doy_median_frac.get('doy_earlySummer')), strategy=segment_find_strategy),
        extrapolate_max_days
    ).regexpRename('$', '_earlySummer').regexpRename('^', 'spectral_')

    peak_summer_fit = fun_refl.fitT(
        seasonal_doy.select('doy_peakSummer'), year,
        segments.find_by_date(ee.Date(doy_median_frac.get('doy_peakSummer')), strategy=segment_find_strategy),
        extrapolate_max_days
    ).regexpRename('$', '_peakSummer').regexpRename('^', 'spectral_')

    late_summer_fit = fun_refl.fitT(
        seasonal_doy.select('doy_lateSummer'), year,
        segments.find_by_date(ee.Date(doy_median_frac.get('doy_lateSummer')), strategy=segment_find_strategy),
        extrapolate_max_days
    ).regexpRename('$', '_lateSummer').regexpRename('^', 'spectral_')

    end_snowfree_fit = fun_refl.fitT(
        seasonal_doy.select('doy_endSnowfree'), year,
        segments.find_by_date(ee.Date(doy_median_frac.get('doy_endSnowfree')), strategy=segment_find_strategy),
        extrapolate_max_days
    ).regexpRename('$', '_endSnowfree').regexpRename('^', 'spectral_')

    return (start_snowfree_fit.addBands(early_summer_fit)
                              .addBands(peak_summer_fit)
                              .addBands(late_summer_fit)
                              .addBands(end_snowfree_fit)
                              .int16()
                              .set('system:index', str(year)))


# 2.1 ----- TOPOGRAPHIC CORRECTION -----

def precompute_topo_c(refl_by_year, slp, ic, seasons, bands_topo_corr, roi,
                                crs, scale, tile_scale):
    """
    Precompute SCSc topographic correction c values for all years × seasons × bands in a
    single batched getInfo() call.

    This function batches all year × season × band linearFit computations into a single 
    ee.List and materializes them with ONE getInfo() call. GEE evaluates all independent 
    reduceRegions calls in parallel server-side — estimated ~15–60 s for 
    40 × 35 = 1,400 fits (vs. 40 × 3–8 s = 2–5 min for serial per-year calls).

    The returned c_values dict maps (year, season, band_name) → Python float. Pass it to
    apply_topo_correction_from_c() to apply the correction per year with no additional
    getInfo() calls and no lazy linearFit chains in the export graph.

    NOTE: The single getInfo() request serializes all N_fits computation graph nodes. For
    40 years × 35 fits = 1,400 nodes this can be a large payload (~MB range). If GEE raises
    a serialization error, process a subset of years and call this function multiple times.

    Parameters
    ----------
    refl_by_year     : dict {year (int) : ee.Image} — seasonal reflectance per year, from
                       compute_seasonal_reflectance(). Keys must be Python ints.
    slp              : ee.Image — Copernicus 30m slope (from load_static_predictors)
    ic               : ee.Image — illumination condition (from load_ccdc_data)
    seasons          : list[str] — season names, e.g. ['startSnowfree', ..., 'endSnowfree']
    bands_topo_corr  : list[str] — band prefixes, e.g. ['spectral_blue_', ...]
    roi              : ee.FeatureCollection — compact feature(s) for the linearFit
    crs, scale       : projection parameters
    tile_scale       : int — GEE tile scale

    Returns
    -------
    dict : {(year, season, band_name) : float}
        band_name is the full band name, e.g. 'spectral_blue_startSnowfree'
    """
    fit_list = []    # lazy ee.Dictionary objects, one per (year, season, band)
    fit_keys = []    # (year, season, band_name) in the same order as fit_list

    for year in sorted(refl_by_year.keys()):
        refl = refl_by_year[year]

        for season in seasons:

            ic_season   = ic.select('IC_'   + season).divide(10000)
            cosZ_season = ic.select('cosZ_' + season).divide(10000)
            cosS_season = ic.select('cosS_' + season).divide(10000)

            refl_ic = ee.Image(refl.divide(10000)
                                   .addBands(ic_season)
                                   .addBands(cosZ_season.rename('cosZ'))
                                   .addBands(cosS_season.rename('cosS'))
                                   .addBands(slp.rename('slope')))
            refl_ic = ee.Image(refl_ic.updateMask(refl_ic.select('slope').gte(5)))
            refl_ic = ee.Image(refl_ic.updateMask(refl_ic.select('IC_' + season).gt(0)))

            for band_prefix in bands_topo_corr:
                band = band_prefix + season
                fit_fc = refl_ic.select(ee.List(['IC_' + season, band])).reduceRegions(
                    collection=roi,
                    reducer=ee.Reducer.linearFit(),
                    crs=crs,
                    scale=scale,
                    tileScale=tile_scale
                )
                fit_list.append(ee.Feature(fit_fc.first()).toDictionary())
                fit_keys.append((year, season, band))

    # Single getInfo() materializes all year × season × band fits.
    # GEE evaluates the independent reduceRegions calls in parallel server-side.
    all_fits = ee.List(fit_list).getInfo()

    # Parse c = offset / scale for each (year, season, band).
    # Degenerate cases same as apply_topo_correction_feat_eecu:
    #   (1) No valid pixels → fit dict is empty → c = 0.0 → fully-masked SCSc output →
    #       unmask(refl.select(band), False) restores original uncorrected reflectance. ✓
    #   (2) scale == 0 (zero IC variance) → c = 0.0 → Lambertian correction. ✓
    c_values = {}
    for (year, season, band), fit in zip(fit_keys, all_fits):
        offset    = (fit or {}).get('offset') # If fit is None, return empty diectionary, .get() returns None
        scale_val = (fit or {}).get('scale') # If fit is None, return empty diectionary, .get() returns None
        if offset is None or scale_val is None or scale_val == 0:
            c_values[(year, season, band)] = 0.0
        else:
            c_values[(year, season, band)] = offset / scale_val

    return c_values


def apply_topo_correction_from_c(refl, slp, ic, seasons, bands_topo_corr, year, c_values):
    """
    Apply SCSc topographic correction using pre-computed c values.

    For each season × band, the c value is embedded as a constant ee.Number — no lazy
    reduceRegions chain remains in the export graph, keeping per-tile EECU at 1–10.

    No crs/scale/tile_scale parameters are needed here — all reduceRegions work was done
    in precompute_topo_c; this function only builds lazy expression graph nodes.

    Parameters
    ----------
    refl             : ee.Image — seasonal reflectance for 'year'
    slp              : ee.Image — Copernicus 30m slope
    ic               : ee.Image — illumination condition
    seasons          : list[str]
    bands_topo_corr  : list[str] — band prefixes
    year             : int — key used to look up c_values[(year, season, band)]
    c_values         : dict — returned by precompute_topo_c()

    Returns
    -------
    ee.Image — lazy multiband image with same band names as input 'refl'
    """
    refl_topo_corr = []

    for season in seasons:

        ic_season   = ic.select('IC_'   + season).divide(10000)
        cosZ_season = ic.select('cosZ_' + season).divide(10000)
        cosS_season = ic.select('cosS_' + season).divide(10000)

        refl_ic = ee.Image(refl.divide(10000)
                               .addBands(ic_season)
                               .addBands(cosZ_season.rename('cosZ'))
                               .addBands(cosS_season.rename('cosS'))
                               .addBands(slp.rename('slope')))
        refl_ic = ee.Image(refl_ic.updateMask(refl_ic.select('slope').gte(5)))
        refl_ic = ee.Image(refl_ic.updateMask(refl_ic.select('IC_' + season).gt(0)))

        season_bands = [b + season for b in bands_topo_corr]
        refl_SCSccorr = []

        for band in season_bands:

            c = ee.Number(c_values[(year, season, band)])

            SCSc_output = refl_ic.expression(
                "((image * (cosB * cosZ + cvalue)) / (ic + cvalue))", {
                    'image':  refl_ic.select(ee.List([band])),
                    'ic':     refl_ic.select('IC_' + season),
                    'cosB':   refl_ic.select('cosS'),
                    'cosZ':   refl_ic.select('cosZ'),
                    'cvalue': c
                }
            )
            SCSc_output = SCSc_output.multiply(10000).unmask(refl.select(band), False)
            refl_SCSccorr.append(ee.Image(SCSc_output))

        refl_SCSccorr = (ee.ImageCollection.fromImages(refl_SCSccorr)
                         .toBands()
                         .regexpRename('^[0-9]_', ''))
        refl_topo_corr.append(refl_SCSccorr)

    return (ee.ImageCollection.fromImages(refl_topo_corr)
            .toBands()
            .int16()
            .regexpRename('^[0-9]_', ''))


# ==========================
# 3. MODEL PREDICTORS ======
# ==========================

# 3.0 ----- AUXILIARY PREDICTORS -----

def compute_aux_predictors(modeled_reflectance, segments, year, crs, scale,
                           texture_bands, texture_metrics, texture_radius,
                           segment_find_strategy, extrapolate_max_days):
    """
    Compute NDVI texture and NDVI slope auxiliary predictors.

    Port of 06_2 logic. Texture is computed on topo-corrected reflectance + normalized indices.
    NDVI slope is derived from the CCDC segment active at mid-summer.
    Both outputs have gap-filled by focal mean and then blended with the original values.

    Parameters
    ----------
    refl, refl_tc        : ee.Image — raw seasonal reflectance and topo-corrected version
    segments             : Segments — from load_ccdc_data(); used to get NDVI slope
    year                 : int — used to define the mid-summer date for segment selection
    crs, scale           : projection parameters
    texture_bands        : list[str] — band(s) to compute texture on, e.g. ['spectral_NDVI_peakSummer']
    texture_metrics      : list[str] — texture metric regex, e.g. ['.*_var'] (variance)
    texture_radius       : int — neighbourhood radius in pixels (1 = 3x3 window)
    segment_find_strategy: str — 'next' or 'closest'
    extrapolate_max_days : int — days to extrapolate beyond CCDC segment endpoints for NDVIslope

    Returns
    -------
    dict with keys:
      'texture' : ee.Image — NDVI texture (uint16), gap-filled and blended
      'slope'   : ee.Image — NDVI slope from CCDC segment (int32), gap-filled and blended
    """

    # NDVI texture (variance in neighbourhood window determined by texture_radius)
    texture = fun_refl.calculateTexture(
        modeled_reflectance, texture_radius, texture_bands, texture_metrics, crs, scale
    )
    texture_fill = texture.focalMean(
        radius=scale * 3 // 2, kernelType='square', units='meters', iterations=1
    )
    texture = texture_fill.blend(texture).uint16()

    # NDVI slope from CCDC segment at mid-summer
    date = ee.Date.fromYMD(year, 7, 31)
    segment = segments.find_by_date(date, strategy=segment_find_strategy)
    slope = fun_refl.NDVIslope(segment, extrapolate_max_days)
    slope_fill = slope.focalMean(
        radius=scale * 3 // 2, kernelType='square', units='meters', iterations=1
    )
    slope = slope_fill.blend(slope).int32().setDefaultProjection(crs=crs, scale=scale)

    return {'texture': texture, 'slope': slope}


# 3.1 ----- SPECTRAL PREDICTORS -----

def compute_spectral_predictors(modeled_reflectance, phenology):
    """
    Compute annual spectral summary predictors (mean, median, range, seasonal change rates).

    Port of 08_1/08_2 spectral predictor block. Calls calcAnnualSpectral for 11 bands/indices.

    Parameters
    ----------
    modeled_reflectance : ee.Image — topo-corrected reflectance + normalized indices
                          (refl_tc.addBands(refl.select('.*_ND.*|.*NBR.*')))
    phenology           : ee.Image — seasonal DOY mosaic (year-independent)

    Returns
    -------
    ee.Image — lazy multiband image with annual spectral summaries per band/index
    """
    # Per-band seasonal summaries (mean, median, range); change rates computed below
    summaries = (fun_misc.calcAnnualSpectral(modeled_reflectance, 'blue')
                 .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'green'))
                 .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'red'))
                 .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'NIR'))
                 .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'SWIR1'))
                 .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'SWIR2'))
                 .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'EVI2b'))
                 .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'NBR'))
                 .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'NDMI'))
                 .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'NDVI'))
                 .addBands(fun_misc.calcAnnualSpectral(modeled_reflectance, 'NDWI')))

    # Change rates: all 11 bands × 4 season pairs in 4 batched operations (vs. 44 individual nodes)
    change_rates = fun_misc.calcAllChangeRates(modeled_reflectance, phenology)

    return summaries.addBands(change_rates)


# 3.2 ----- TREE COVER -----

def get_tree_cover(year):
    """
    Return year-appropriate tree cover image.

    Uses Python if/else (year is always a Python int, so no server-side conditional needed).
    For 2020 uses a pre-processed asset; for all other years uses Hansen 2000 treecover.

    Returns
    -------
    ee.Image — lazy tree cover image
    """
    if year == 2020:
        return ee.Image("projects/arctic-biomass-mapping/assets/predictors/tree_cover_2020")
    return (ee.Image("UMD/hansen/global_forest_change_2022_v1_10")
            .select(['treecover2000'], ['cover'])
            .regexpRename('^', 'trees_')
            .addBands(
                ee.Image("UMD/hansen/global_forest_change_2022_v1_10")
                  .select('treecover2000')
                  .gt(0)
                  .rename('trees_presence')
            )
            .unmask(0))


# 3.3 ----- TREE COVER -----

def assemble_predictors(spectral, static, tree_cover, texture, crs, scale):
    """
    Assemble the full predictor image for model classification.

    Call this ONCE before the MC loop and pass the result to run_binary_mc_ic() /
    run_continuous_mc_ic(). Assembling inside the MC loop would replicate the entire
    predictor graph 100× per IC, causing 'Encoded object too large' errors when
    sampleRegions() or reduceRegions() serializes the combined computation graph.

    All possible predictor bands are always included — topographic, ecoregions, and all
    others. GEE classifiers select features by name from the input image and silently
    ignore extra bands, so a single predictor image is safe for both binary models
    (trained without ecoregions) and continuous models (trained with ecoregions).

    Parameters
    ----------
    spectral   : ee.Image — from compute_spectral_predictors()
    static     : dict — from load_static_predictors()
    tree_cover : ee.Image — from get_tree_cover()
    texture    : ee.Image — from compute_aux_predictors()['texture']
    crs, scale : projection parameters

    Returns
    -------
    ee.Image — lazy predictor stack ready to classify
    """
    lat_lon = ee.Image.pixelLonLat()

    return (spectral
        .addBands(lat_lon)
        .addBands(static['permafrost'])
        .addBands(texture)
        .addBands(static['wte'].reduceResolution(reducer=ee.Reducer.mode(), maxPixels=1000))
        .addBands(ee.Image.constant(0).rename('world_terrestrial_ecosystems_other')
                    .setDefaultProjection(crs=crs, scale=scale)
                    .reduceResolution(reducer=ee.Reducer.mode(), maxPixels=1000))
        .addBands(tree_cover)
        .addBands(static['zones'].reduceResolution(reducer=ee.Reducer.mode(), maxPixels=1000))
        .addBands(static['topographic'])
        .addBands(static['ecoregions'].reduceResolution(reducer=ee.Reducer.mode(), maxPixels=1000)))


# ================================
# 4. MONTE CARLO ITERATIONS ======
# ================================
# Run with Monte Carlo iterations as the outer loop, and years as the inner loop
#
# WHY MC-OUTER: In the year-outer approach, run_continuous / run_binary is called once per year,
# producing 100 blob references (MC models) per year. With 40 years stacked,
# the export graph holds 100 × 40 = 4,000 blob references. GEE inlines blob string content
# per tile, causing OOM for continuous trees (~2.55 MB × 4,000 ≈ 10 GB per tile).
#
# In the MC-outer approach, a single ee.Blob / ee.Classifier Python object is created once
# per MC and shared across all year classify() calls. GEE's graph deduplicates shared Python
# objects, so the export graph holds exactly 100 blob references regardless of year count.

def run_continuous_mc(ds_type, model_version, mc_list,
                                predictors_by_year, year_list, tree_dir):
    """
    Build a 100-image MC IC for continuous predictions across all years (MC-outer).

    Parameters
    ----------
    ds_type            : str — 'total' or 'woody'
    model_version      : str — e.g. 'v20240514'
    mc_list            : list[int] — MC iteration numbers, e.g. list(range(1, 101))
    predictors_by_year : dict[int, ee.Image] — pre-built, short-renamed predictor image per year
    year_list          : list[int] — ordered list of years to process
    tree_dir           : str — GCS subdirectory for tree files

    Returns
    -------
    ee.ImageCollection — 100 images; each has year-banded uint16 predictions with bands
                         named {ds_type}_continuous_{year}
    """
    mc_images = []
    for mc in mc_list:
        path = ('gs://arctic_biomass_mapping_models/' + model_version
                + '/mc/continuous/' + tree_dir + '/'
                + ds_type + '_continuous_formatted_gee_short_'
                + model_version + '_' + str(int(mc)) + '.txt')
        # One blob + one classifier per MC; shared across all year classify() calls below
        mod = ee.Classifier.decisionTreeEnsemble(
            ee.List([ee.Blob(path).string()])
        ).setOutputMode('REGRESSION')

        year_bands = []
        for year in year_list:
            pred = (predictors_by_year[year]
                    .classify(mod, 'predicted')
                    .round().uint16()
                    .rename(ds_type + '_continuous_' + str(year)))
            year_bands.append(pred)

        mc_img = (ee.ImageCollection(year_bands)
                  .toBands()
                  .regexpRename('^[0-9]+_', ''))
        mc_images.append(mc_img)

    return ee.ImageCollection(mc_images)


def run_binary_mc(ds_type, model_version, mc_list,
                            predictors_by_year, year_list, tree_dir):
    """
    Build a 100-image MC IC of raw presence probability for all years (MC-outer).

    Returns raw prob_presence × 100 (Int8) per year-band per MC. The same IC can be
    passed to both aggregate_binary_mc and aggregate_probability_mc,
    avoiding running the model twice.

    Parameters
    ----------
    ds_type            : str — 'total' or 'woody'
    model_version      : str — e.g. 'v20240514'
    mc_list            : list[int]
    predictors_by_year : dict[int, ee.Image] — pre-built, short-renamed predictor image per year
    year_list          : list[int]
    tree_dir           : str — GCS subdirectory for tree files

    Returns
    -------
    ee.ImageCollection — 100 images; each has year-banded Int8 presence probability × 100
                         with bands named {ds_type}_binary_{year};
                         system:index set to MC number string (e.g. '1' … '100') for
                         per-MC threshold lookup in aggregate_binary_mc
    """
    mc_images = []
    for mc in mc_list:
        path = ('gs://arctic_biomass_mapping_models/' + model_version
                + '/mc/binary/' + tree_dir + '/'
                + ds_type + '_binary_formatted_gee_short_'
                + model_version + '_' + str(int(mc)) + '.txt')
        mod = ee.Classifier.decisionTreeEnsemble(
            ee.List([ee.Blob(path).string()])
        ).setOutputMode('RAW')

        year_bands = []
        for year in year_list:
            results = predictors_by_year[year].classify(mod, 'raw_probability')
            prob = (results.arrayReduce(ee.Reducer.mean(), ee.List([0]))
                          .arrayGet(0)
                          .multiply(100).toInt8()
                          .rename(ds_type + '_binary_' + str(year)))
            year_bands.append(prob)

        mc_img = (ee.ImageCollection(year_bands)
                  .toBands()
                  .regexpRename('^[0-9]+_', '')
                  .set('system:index', str(int(mc))))
        mc_images.append(mc_img)

    return ee.ImageCollection(mc_images)


# =========================
# 5. MAP AGGREGATION ======
# =========================

def aggregate_continuous_mc(mc_ic, ds_type, year_list):
    """
    Reduce a continuous MC IC to a 3-band-per-year percentile summary.

    Parameters
    ----------
    mc_ic     : ee.ImageCollection — from run_continuous_mc(); 100 images,
                each with year-banded predictions named {ds_type}_continuous_{year}
    ds_type   : str
    year_list : list[int]

    Returns
    -------
    ee.Image — 3 bands per year, uint16:
               {ds_type}_biomass_{year}      (p50 across MC)
               {ds_type}_biomass_lwr_{year}  (mean of p2, p3)
               {ds_type}_biomass_upr_{year}  (mean of p97, p98)
    """
    cont_pct = mc_ic.reduce(ee.Reducer.percentile([2, 3, 50, 97, 98]))
    cont_pct = cont_pct.updateMask(cont_pct.mask().gt(0))

    year_outputs = []
    for year in year_list:
        band_base = ds_type + '_continuous_' + str(year)
        median = cont_pct.select(band_base + '_p50').rename(ds_type + '_biomass_' + str(year))
        lwr = ee.ImageCollection([
            cont_pct.select(band_base + '_p2').rename(ds_type + '_biomass_lwr_' + str(year)),
            cont_pct.select(band_base + '_p3').rename(ds_type + '_biomass_lwr_' + str(year))
        ]).mean()
        upr = ee.ImageCollection([
            cont_pct.select(band_base + '_p97').rename(ds_type + '_biomass_upr_' + str(year)),
            cont_pct.select(band_base + '_p98').rename(ds_type + '_biomass_upr_' + str(year))
        ]).mean()
        year_outputs.append(median.addBands(lwr).addBands(upr))

    return (ee.ImageCollection(year_outputs)
            .toBands()
            .regexpRename('^[0-9]+_', '')
            .round().uint16())


def aggregate_binary_mc(mc_ic, thresholds):
    """
    Apply per-MC thresholds to a binary MC IC and take mode across MC iterations.

    Applies each MC's probability threshold to all year-bands of that MC's image
    simultaneously (gte() with a scalar threshold broadcasts to all bands), then takes
    the mode across 100 MC images per year-band.

    Parameters
    ----------
    mc_ic      : ee.ImageCollection — from run_binary_mc(); each image has
                 year-banded prob_presence × 100 (Int8) and system:index = MC number string
    thresholds : ee.FeatureCollection — per-MC threshold; must have numeric 'mc' property
                 and 'threshold' property in [0, 1]

    Returns
    -------
    ee.Image — 1 band per year, byte:
               {ds_type}_presence_absence_{year}
    """
    def apply_threshold(mc_image):
        mc = ee.String(mc_image.get('system:index'))
        threshold_val = (ee.Number(
            thresholds.filter(ee.Filter.eq('mc', ee.Number.parse(mc))).first().get('threshold')
        ).multiply(100).round())
        return mc_image.gte(threshold_val).byte()

    mode_img = mc_ic.map(apply_threshold).reduce(ee.Reducer.mode())
    mode_img = mode_img.updateMask(mode_img.mask().gt(0))
    return (mode_img
            .regexpRename('_mode$', '')
            .regexpRename('_binary_', '_presence_absence_')
            .byte())


def aggregate_probability_mc(mc_ic, ds_type, year_list):
    """
    Reduce a binary MC IC to a 3-band-per-year probability percentile summary.

    Pass the same ee.ImageCollection returned by run_binary_mc() that was used
    for aggregate_binary_mc — no rerunning of the model is needed.

    Parameters
    ----------
    mc_ic     : ee.ImageCollection — from run_binary_mc(); each image has
                year-banded prob_presence × 100 (Int8)
    ds_type   : str
    year_list : list[int]

    Returns
    -------
    ee.Image — 3 bands per year, byte:
               {ds_type}_probability_{year}      (p50 × 100)
               {ds_type}_probability_lwr_{year}  (mean of p2, p3 × 100)
               {ds_type}_probability_upr_{year}  (mean of p97, p98 × 100)
    """
    prob_pct = mc_ic.reduce(ee.Reducer.percentile([2, 3, 50, 97, 98]))
    prob_pct = prob_pct.updateMask(prob_pct.mask().gt(0))

    year_outputs = []
    for year in year_list:
        band_base = ds_type + '_binary_' + str(year)
        median = prob_pct.select(band_base + '_p50').rename(ds_type + '_probability_' + str(year))
        lwr = ee.ImageCollection([
            prob_pct.select(band_base + '_p2').rename(ds_type + '_probability_lwr_' + str(year)),
            prob_pct.select(band_base + '_p3').rename(ds_type + '_probability_lwr_' + str(year))
        ]).mean()
        upr = ee.ImageCollection([
            prob_pct.select(band_base + '_p97').rename(ds_type + '_probability_upr_' + str(year)),
            prob_pct.select(band_base + '_p98').rename(ds_type + '_probability_upr_' + str(year))
        ]).mean()
        year_outputs.append(median.addBands(lwr).addBands(upr))

    return (ee.ImageCollection(year_outputs)
            .toBands()
            .regexpRename('^[0-9]+_', '')
            .round().byte())