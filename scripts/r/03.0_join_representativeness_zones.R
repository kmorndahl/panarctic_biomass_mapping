################################################################################
################################################################################

# DESCRIPTION:

# This script:
#   - summarizes plot sizes (mean and mode) across PFTs for aggregate datasets (e.g. total, nontree, woody)
#   - combines synthesis dataset and GEE representativeness analysis
#   - extracts vector-based predictors: arctic zone, ecoregion

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

# NOTES:
# - For water/snow/wetland fraction, number of pixels is higher than expected
#   - Expectation:
#     - 30x30m = 3x3 10m sentinel pixels = 9 total pixels
#     - 90x90m = 9x9 10m sentinel pixels = 81 total pixels
# - Different from expected because WorldCover data is in different projection with rectangular pixels and slightly smaller pixel resolution

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

rm(list=ls())
require(data.table)
require(tidyverse)
library(dtplyr)
require(sf)
require(R.utils)
require(tibble)

# 1.2 Parameters ---------------------------------------------------------------

rep_version = 'v20231110'
out_version = 'v20240215'

plot_CRS = 4326

synth_dir = 'data/01_synthesis_dataset/'
rep_dir = paste0('output/03_representativeness_analysis/csv_files/')
nonveg_dir = 'data/01_synthesis_dataset/'
zone_dir = 'data/000_gis/'
out_dir_joined = 'output/03_representativeness_analysis/joined/'

ds_name_list = c('plant', 'woody_plant')

# 1.3 Functions ----------------------------------------------------------------

find_mode = function(x) {
  u = unique(na.omit(x))
  tab = tabulate(match(x, u))
  u[tab == max(tab)][1] # If there are two modes, grab the first one
}

# 1.4 Read in and tidy data ----------------------------------------------------

for(ds_name in ds_name_list){

  # Read in dataset
  if(ds_name == 'all'){
    ds.synth = fread(paste0(synth_dir, 'the_arctic_plant_aboveground_biomass_synthesis_dataset_v1.2.csv'))
  }else{
    ds.synth = fread(paste0(synth_dir, 'arctic_tundra_total_', ds_name, '_biomass_synthesis_dataset.csv'))
  }
  
  # Tidy dataset names
  ds_name = gsub('_plant', '', ds_name)
  if(ds_name == 'plant') ds_name = 'total'

  # Read in full synthesis dataset (necessary for plot size information)
  ds.synth.all = fread(paste0(synth_dir, 'the_arctic_plant_aboveground_biomass_synthesis_dataset_v1.2.csv'))

  # Rename if necessary
  if(ds_name %in% c('plant', 'nontree', 'woody')){
    ds.synth = ds.synth %>% rename("biomass_density_gm2" = "plant_biomass_density_gm2")
  }

  # Remove unnecessary columns
  ds.synth = ds.synth %>% select(-any_of(c('gdd_degC', 'mat_degC', 'map_mm', 'community_measured', 'plants_measured', 'nontree_plants_measured', 'vascular_measured', 'woody_measured', 'biomass_dry_weight_g', 'method', 'vegetation_description', 'site_description', 'notes', 'bioclim_zone')))
  if(ds_name != 'all') ds.synth = ds.synth %>% select(-any_of(c('pft')))
  
  # Save original data.table size
  ds.orig.size = nrow(ds.synth)

  # Read in representativeness files
  # Using .csv files because complicated issue reading in .kml files... https://github.com/r-spatial/sf/issues/499
  rep.files = list.files(rep_dir, full.names = T)
  ds.rep = do.call("rbind", lapply(rep.files, fread))
  ds.rep = ds.rep %>% rename_at(vars(NBR_mean:wetland_pixel_count), ~ paste0("rep_", .x)) %>% as.data.table() # Label representativeness columns

  # Read in nonvegetated points
  sf.nonveg = st_read(paste0(nonveg_dir, 'barren_built_snow_ice_water_random_pts_scale30m.shp'))
  sf.nonveg = st_transform(sf.nonveg, plot_CRS)

  # Read in zones/ecoregions
  arctic_zones_3 = sf::st_read(paste0(zone_dir, 'arctic_three_lat_zones_laea.shp'))
  arctic_ecoregions = sf::st_read(paste0(zone_dir, 'ecoregions_arctic.shp'))

  # ==============================================================================
  # 2. SYNTHESIS DATA ============================================================
  # ==============================================================================

  # 2.1 Add plot size to datasets ------------------------------------------------

  if(ds_name %in% c('plant', 'nontree', 'woody')){
    ds.synth = merge(ds.synth, ds.synth.all[, c("plot_code", "plot_area_m2", "pft")], by = "plot_code", all.x = TRUE)
  }

  # Some some plots have varying plot sizes depending on PFT, calculate mean and mode

  if(ds_name == 'all'){
    ds.synth = ds.synth %>%
      group_by(plot_code) %>%
      mutate(
        plot_area_m2_mode = find_mode(plot_area_m2),
        plot_area_m2_mean = mean(plot_area_m2, na.rm = TRUE) # For 'all' dataset, ignore NAs - these are PFTs that were not measured
      ) %>%
      as.data.table()
  }else if(ds_name %in% c('plant', 'nontree', 'woody')){
    ds.synth = ds.synth[ds.synth$pft != 'lichen',] %>%
      group_by(plot_code) %>%
      summarise(
        across(-plot_area_m2, first),
        plot_area_m2_mode = find_mode(plot_area_m2),
        plot_area_m2_mean = mean(plot_area_m2)
      ) %>%
      dplyr::select(-any_of(c('pft'))) %>%
      as.data.table()
  }else{
    ds.synth = ds.synth %>%
      group_by(plot_code) %>%
      summarise(
        across(-plot_area_m2, first),
        plot_area_m2_mode = find_mode(plot_area_m2),
        plot_area_m2_mean = mean(plot_area_m2)
      ) %>%
      dplyr::select(-any_of(c('pft'))) %>%
      as.data.table()
  }

  # Check
  if(nrow(ds.synth) != ds.orig.size){stop('Synthesis data.table different length after adding plot area, check summarise function')}

  # 2.1 Add non-vegetated observations to synthesis dataset ----------------------

  # Get non-vegetated observations only
  nonveg = ds.rep[ds.rep$dataset_id == 'nonveg',]
  nonveg.size = nrow(nonveg)

  # Remove representativeness columns
  cols = grep("NBR|NDMI|NDVI|water|snow|wetland", names(nonveg), ignore.case=TRUE, value = TRUE)
  nonveg = nonveg[, !..cols]

  # Get plot ID
  nonveg = nonveg[, c("site_id", "plot_id") := tstrsplit(plot_code, "_", fixed = TRUE)] # Split 'plot_code' column into 'site_id' and 'plot_id' columns by '_'
  nonveg$plot_id = as.numeric(nonveg$plot_id)

  # Join lat/long data
  sf.nonveg = cbind(sf.nonveg %>% st_drop_geometry(), # Drop geometry
                    data.frame(st_coordinates(sf.nonveg))) %>% # Get geometry as coordinates
    dplyr::select(-c(Map, acceptable)) # Remove unneeded columns
  nonveg = merge(nonveg, sf.nonveg, all.x = TRUE, by.x="plot_id", by.y="id")
  nonveg = nonveg %>% rename("latitude" = "Y", "longitude" = "X") %>% as.data.table()

  # Populate nonveg fields to match synthesis plot data
  nonveg$contributor = "nonveg"
  nonveg$country = "nonveg"
  nonveg$locale = "nonveg"
  nonveg$year = 2021
  nonveg$month = 7
  nonveg$day = 31
  nonveg$citation = "nonveg"
  nonveg$biomass_density_gm2 = 0
  nonveg$plot_area_m2_mode = 900
  nonveg$plot_area_m2_mean = 900
  if(ds_name == 'all'){
    nonveg$plot_area_m2 = 900
    nonveg$pft = "total"
  }

  # Add to synthesis dataset
  ds.synth = rbind(ds.synth, nonveg)

  # Size after adding nonveg
  ds.size = nrow(ds.synth)

  # Check
  if(ds.size != (ds.orig.size + nonveg.size)){stop('After adding nonveg, synthesis data.table size does not match original size + nonveg size, check add nonveg code')}

  # ==============================================================================
  # 3. COMBINE DATA ==============================================================
  # ==============================================================================

  # Join synthesis data with representativeness datasets
  ds = merge(ds.synth, ds.rep, by = c("dataset_id", "citation_short", "site_code","plot_code", "coord_type"), all.x = TRUE)

  # Check
  if(nrow(ds) != ds.size){stop('Synthesis data.table size incorrect after joining, check join code')}

  # ==============================================================================
  # 4. EXTRACT ZONES/ECOREGIONS ==================================================
  # ==============================================================================

  # 4.1 Tidy spatial data ----------------------

  # Convert to spatial
  sf = st_as_sf(as.data.frame(ds) %>% drop_na(latitude), coords = c("longitude", "latitude"), crs = plot_CRS)

  # Reproject - need to reproject plots instead of arctic zones because reprojecting arctic zones produces invalid polygons
  sf = sf::st_transform(sf, sf::st_crs(arctic_zones_3))
  arctic_ecoregions = sf::st_transform(arctic_ecoregions, sf::st_crs(arctic_zones_3))

  # Subset zone datasets
  arctic_zones_3 = arctic_zones_3 %>% dplyr::select(FIRST_zone)
  arctic_ecoregions = arctic_ecoregions %>% dplyr::select(ECO_NAME)

  # 4.2 Extract zones ----------------------

  # Some points fall just outside of arctic zones, apply within distance filter
  # Those outside of the distance filter are assigned 'Subarctic'
  sf = st_join(sf, arctic_zones_3, join = st_is_within_distance, dist = 500, left = TRUE)

  # Remove duplicates
  # Can have duplicates if more than one zone polygon within 500m
  # st_join data is sorted alphabetically so when choosing which observation to retain from duplicates it will prioritize in this order:
  #   High Arctic; Low Arctic; Oro Arctic, Sub Arctic
  if(ds_name == 'all') sf = sf[!duplicated(sf[c("plot_code","pft")]),] else sf = sf[!duplicated(sf$plot_code),]

  # Assign 'subarctic' class wherever plot did not intersect one of the arctic zones
  sf$FIRST_zone[is.na(sf$FIRST_zone)] = 'Sub Arctic'

  # 4.3 Extract ecoregions ----------------------

  # Some points fall outside of ecoregions, because ecoregions cover our whole study area we can apply a nearest feature join
  sf = st_join(sf, arctic_ecoregions, join = st_nearest_feature, left = TRUE)

  # 4.4 Tidy ----------------------

  # Rename zone columns
  sf = sf %>% rename(zone_3 = FIRST_zone, ecoregion = ECO_NAME)

  # Convert back to WGS84
  sf = sf::st_transform(sf, plot_CRS)

  # Return to dataframe, preserving lat/long information
  ds.final = cbind(st_drop_geometry(sf), st_coordinates(sf)) %>% rename(latitude = Y, longitude = X)

  # Check
  if(nrow(ds.final) != nrow(ds)){stop('Synthesis data.table different length after joining zones and ecoregions, check st_join')}

  # ==============================================================================
  # 5. SAVE ======================================================================
  # ==============================================================================

  # Check NAs
  # NAs for representativeness data are okay -- they will be filtered out later
  # NAs for day are okay -- the will be filled in later
  # For the full synthesis dataset ('all') NAs for biomass and plot area are okay -- these are unmeasured PFTs
  colSums(is.na(ds.final))

  # Write out
  fwrite(ds.final, paste0(out_dir_joined, 'arctic_tundra_', ds_name, '_biomass_synthesis_representativeness_joined_', out_version, '.csv'))

}


