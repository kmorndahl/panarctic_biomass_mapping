"""

DESCRIPTION: Generate site representativeness summary statistics

AUTHOR: Melissa Rose, NAU
DATE: 1-25-2023

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

# 1.0 ----- EXISTING MODULES -----

# Import required modules
from utils import field_data as fun_field_data
from utils import sentinel as fun_sentinel

# 1.1 ----- PARAMETERS -----

# Output options
version = 'v20240215'
export_loc = 'drive'  # Location to export the representativeness summary statistics (options: drive, asset)
asset_folder = 'arctic_biomass/'  # If 'asset' is specified above, specify your asset folder (asset name is constructed automatically)
download_folder = 'arctic_biomass_synthesis_representativeness'  # Specify folder in Google drive to download output csvs and shapefiles to
export_properties = [
    'dataset_id', 'citation_short', 'site_code', 'plot_code', 'coord_type',
    'NBR_mean', 'NBR_stdDev',
    'NDMI_mean', 'NDMI_stdDev',
    'NDVI_mean', 'NDVI_stdDev',
    'nbr_cv', 'ndmi_cv', 'ndvi_cv',
    'snow_fraction', 'snow_pixel_count',
    'water_fraction', 'water_pixel_count',
    'wetland_fraction', 'wetland_pixel_count',
    'total_pixel_count'
]  # Set properties to export - must set list of export properties manually otherwise get 'Too many concurrent aggregations' error

# 1.2 ----- DATASET INPUTS -----

# Import MGRS tiles
tiles = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/ROIs/tiles_mgrs_s2_arctic_update').select('epsg')

# Import coordinates for all datasets and filter to obtain a single lat/long value for each plot
field_data = ee.FeatureCollection('projects/arctic-biomass-mapping/assets/field_data/arctic_tundra_biomass_synthesis_plots_' + version).distinct('plot_code')

# Get list of unique dataset IDs for mapping
datasets = field_data.distinct('dataset_id')

def get_dataset_id(f):
    return ee.Feature(f).getString('dataset_id')

datasets = datasets.toList(100).map(get_dataset_id)
datasets = datasets.getInfo()

# Report data results
print('Field data:', field_data.getInfo())
print('All datasets', datasets)

# Select which dataset to process (set to datasets to process all datasets but export each in an individual .csv and .kml)
select_dataset = datasets  # Output individual .csv per dataset

# 1.3 ----- SENTINEL-2 INPUTS -----

# Define start and end dates for Sentinel-2 cloudless composite images
#   (Image collection will be further filtered to select images only during peak greenness)
start_date = '2019-07-01'
end_date = '2021-08-31'

# Define bands for summary statistics
#   Available bands: ['Blue','Green','Red','RE1','RE2','RE3','NIR','SWIR1','SWIR2']
#   Available indices: ['NDVI', 'kNDVI', 'EVI', 'RGVI', 'NDMI', 'NBR', 'CIRE', 'SVVI']
select_bands = ['NDVI', 'NDMI', 'NBR']

# Input a threshold value for the maximum number of snow/ice or open water pixels allowed in each plot
water_pixel_max = 0

# =========================
# 2. FUNCTION =============
# =========================

feature_collection = field_data
start = start_date
end = end_date
bands = select_bands
max = water_pixel_max
export_type = export_loc
asset = asset_folder
download = download_folder

def representativeness(dataset_id):

    # Step 1: Filter feature collection to specified dataset
    AOI = feature_collection.filter(ee.Filter.equals('dataset_id', dataset_id))

    # Step 2: Create square bounding regions around each plot/site
    AOI_buffer = fun_field_data.buffer_func(dataset_id, AOI)

    # Step 3: Create a 2-year cloudless composite image during peak greenness (mid-July to mid-August) across entire study site
    AOI_S2CC = fun_sentinel.quality_S2_CloudlessComposite(dataset_id, AOI_buffer, start, end, export_type, asset)

    # Step 4: Calculate the mean, standard deviation and coefficient of variation for NDVI, NDMI, and NBR within each plot/site
    AOI_CV = fun_field_data.mean_sd_cv(dataset_id, AOI_buffer, AOI_S2CC, bands)

    # Step 5: Calculate the total number and relative percentage of snow/ice, water, and wetland pixels within each plot/site
    AOI_waterMask = fun_field_data.water_mask(dataset_id, AOI_CV, ee.ImageCollection('COPERNICUS/S2_SR'))  # returns 3-band mask image: 'Snow and ice', 'Open water', and 'Herbaceous wetland'
    AOI_final = fun_field_data.water_summary(dataset_id, AOI_CV, AOI_waterMask)  # returns an updated feature collection with water pixel summaries per plot

    # Step 6: Export results

    # Export as .csv
    task = ee.batch.Export.table.toDrive(
        collection=AOI_final,
        description='arctic_tundra_biomass_synthesis_plots_gee_representativeness_' + dataset_id,
        folder=download,
        fileFormat='CSV',
        selectors=export_properties
    )
    task.start()

    return AOI_final

# =========================
# 3. RUN FUNCTION =========
# =========================

# Set dataset_id to user-specified dataset
ds_list = select_dataset
ds_list_end_idx = len(ds_list)  # Get size of dataset
print('Number of datasets to process:', ds_list_end_idx)
print('Dataset IDs to process:', ds_list)

# Run representativeness function for user-specified dataset
data_rep = list(map(representativeness, ds_list))

raise Exception('stop')
