"""

DESCRIPTION: Group MGRS tiles to facilitate analysis

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

data_version = 'v20231110'

# 1.1 ----- READ IN DATA -----

tiles = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/tiles_mgrs_s2').filter(ee.Filter.gte('lat_deg', 45))  # 45 N or higher

# =========================
# 2. GROUPING =============
# =========================

# 2.0 ----- GROUP BY START DOY -----

# Get list of all start DOYs
doy_list = tiles.aggregate_array('start_doy').distinct().getInfo()

# Loop through start DOYs to group and combine tiles
def group_by_doy(doy):
    tiles_doy = tiles.filter(ee.Filter.eq('start_doy', doy))  # Get tiles with start DOY
    tiles_doy_union = ee.Feature(tiles_doy.union().first())  # Union tiles
    return tiles_doy_union.copyProperties(tiles_doy.first(), ['end_doy', 'epsg', 'start_doy', 'strata', 'zone']).set('name', 'sdoy' + str(doy))

tiles_group_doy = ee.FeatureCollection(list(map(group_by_doy, doy_list)))

# 2.1 ----- GROUP BY EPSG -----

# Get list of all EPSG codes
epsg_list = tiles.aggregate_array('epsg').distinct().getInfo()

# Loop through EPSG codes to group and combine tiles
def group_by_epsg(epsg_num):
    tiles_epsg = tiles.filter(ee.Filter.eq('epsg', epsg_num))  # Get tiles with EPSG code
    tiles_epsg_union = ee.Feature(tiles_epsg.union().first())  # Union tiles
    return tiles_epsg_union.copyProperties(tiles_epsg.first(), ['end_doy', 'epsg', 'start_doy', 'strata', 'zone']).set('name', 'epsg' + str(epsg_num))

tiles_group_epsg = ee.FeatureCollection(list(map(group_by_epsg, epsg_list)))

# 2.2 ----- GROUP BY UTM ZONE -----

# Get list of all UTM zones
utm_list = tiles.aggregate_array('zone').distinct().getInfo()

# Loop through UTM zones to group and combine tiles
def group_by_utm(utm_zone):
    tiles_utm = tiles.filter(ee.Filter.eq('zone', utm_zone))  # Get tiles with UTM zone
    tiles_utm_union = ee.Feature(tiles_utm.union().first())  # Union tiles
    return tiles_utm_union.copyProperties(tiles_utm.first(), ['end_doy', 'epsg', 'start_doy', 'strata', 'zone']).set('name', 'utm' + str(utm_zone))

tiles_group_utm = ee.FeatureCollection(list(map(group_by_utm, utm_list)))

# =========================
# 3. EXPORT ===============
# =========================

task = ee.batch.Export.table.toAsset(
    collection=tiles_group_doy,
    description='tiles_mgrs_s2_grouped_start_doy_45n_tidy',
    assetId='projects/arctic-biomass-mapping/assets/ROIs/tiles_mgrs_s2_grouped_start_doy_45n_tidy'
)
task.start()

task = ee.batch.Export.table.toAsset(
    collection=tiles_group_epsg,
    description='tiles_mgrs_s2_grouped_epsg_45n_tidy',
    assetId='projects/arctic-biomass-mapping/assets/ROIs/tiles_mgrs_s2_grouped_epsg_45n_tidy'
)
task.start()

task = ee.batch.Export.table.toAsset(
    collection=tiles_group_utm,
    description='tiles_mgrs_s2_grouped_utm_zone_45n_tidy',
    assetId='projects/arctic-biomass-mapping/assets/ROIs/tiles_mgrs_s2_grouped_utm_zone_45n_tidy'
)
task.start()

raise Exception('stop')
