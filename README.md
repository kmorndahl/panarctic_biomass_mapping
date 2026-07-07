# Pan-Arctic vegetation biomass mapping

## Project summary

The Arctic is warming faster than anywhere else on Earth, placing tundra ecosystems at the forefront of global climate change. Plant biomass is a fundamental ecosystem attribute that is sensitive to changes in climate, closely tied to ecological function, and crucial for constraining ecosystem carbon dynamics. However, the amount, functional composition, and distribution of plant biomass are only coarsely quantified across the Arctic. Therefore, we developed the first moderate resolution (30 m) maps of live aboveground plant biomass (g m<sup>-2</sup>) and woody plant dominance (%) for the Arctic tundra biome, including the mountainous Oro Arctic. 

### Project correspondence:

Kathleen M. Orndahl
kathleen.orndahl@nau.edu

### Companion manuscript:

If you use this code or data, please cite the companion manuscript:

> Orndahl, K. M., Berner, L. T., Macander, M. J., Arndal, M. F., Alexander, H. D., Humphreys, E. R., ... & Goetz, S. J. (2025). Next generation Arctic vegetation maps: Aboveground plant biomass and woody dominance mapped at 30 m resolution across the tundra biome. *Remote Sensing of Environment*, *323*, 114717. https://doi.org/10.1016/j.rse.2025.114717

```bibtex
@article{orndahl2025arctic,
  author  = {Orndahl, Kathleen M. and Berner, Logan T. and Macander, Matthew J. and Arndal, Marie F. and Alexander, Heather D. and Humphreys, Elyn R. and Loranty, Michael M. and Ludwig, Sarah M. and Nyman, Johanna and Juutinen, Sari and Aurela, Mika and Mikola, Juha and Mack, Michelle C. and Rose, Melissa and Vankoughnett, Mathew R. and Iversen, Colleen M. and Kumar, Jitendra and Salmon, Verity G. and Yang, Dedi and Grogan, Paul and Danby, Ryan K. and Scott, Neal A. and Olofsson, Johan and Siewert, Matthias B. and Deschamps, Lucas and Maire, Vincent and Lévesque, Esther and Gauthier, Gilles and Boudreau, Stéphane and Gaspard, Anna and Bret-Harte, M. Syndonia and Raynolds, Martha K. and Walker, Donald A. and Michelsen, Anders and Kumpula, Timo and Villoslada, Miguel and Ylänne, Henni and Luoto, Miska and Virtanen, Tarmo and Greaves, Heather E. and Forbes, Bruce C. and Heim, Ramona J. and Hölzel, Norbert and Epstein, Howard and Bunn, Andrew G. and Holmes, Robert Max and Natali, Susan M. and Virkkala, Anna-Maria and Goetz, Scott J.},
  title   = {Next generation {Arctic} vegetation maps: Aboveground plant biomass and woody dominance mapped at 30 m resolution across the tundra biome},
  journal = {Remote Sensing of Environment},
  year    = {2025},
  volume  = {323},
  pages   = {114717},
  doi     = {10.1016/j.rse.2025.114717}
}
```

## Repository summary

This repository includes scripts for:
- Aggregating and processing training and validation data
- Aggregating and processing Landsat predictor data, including seasonal reflectance modeling using the Continuous Change Detection and Classification Algorithm
- Modeling plant and woody plant aboveground biomass for the year 2020:
    - Performing a nested cross-validation to assess model accuracy
    - Generating 100 Monte Carlo iterations of final models to assess uncertainty
- Mapping biomass across the Pan-Arctic region:
    - Using Monte Carlo iterations to derive best estimates of and 95% uncertainty intervals for plant and woody plant aboveground biomass
    - Computing woody plant dominance
- Performing a series of analyses assessing the spatial distribution of vegetation biomass across the Arctic, and the influence of topography, disturbance, and climate

Data and output associated with this project have been moved to panarctic_biomass_mapping_data and panarctic_biomass_mapping_output to maintain a lightweight script repository.

Scripts in this repository have been converted to Python from the original JavaScript. To use the original JavaScript code, use release v1.0.

## Repository details

Scripts are prepended with numbers for organization, see details below for details. Some scripts are run on a high performance computing cluster, indicated by "(HPC)"

### R code

000 UTILITIES  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;000.0_yardstick.bias.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;000.1_yardstick_rmpse.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;000.2_reptree_update.R  
00 CREATE REGIONS OF INTEREST AND SUPPLEMENTARY BIOMASS SYNTESIS DATASETS  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;00.0_define_arctic_roi.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;00.1_define_arctic_zones.R  
01 LANDSAT CROSS-SENSOR CALIBRATION  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;01.0_landsat_calibration_export_data.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;01.1_landsat_calibration_clean_data.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;01.2_landsat_calibration_fit_models.R  
02 CREATE DATASETS TO USE IN GEE TO EXTRACT PREDICTOR DATA  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;02.0_gee_data_aggregate_final_plots.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;02.1_gee_data_plot_date_dataset.R  
03 ADD REPRESENTATIVENESS DATA TO BIOMASS DATA  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;03.0_join_representativeness_zones.R  
04 ADD PREDICTOR DATA TO BIOMASS DATA  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;04.0_join_predictors.R  
05 FILTER DATA PLOT LEVEL DATA  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;05.0_filter_plot_level_data.R  
06 AGGREGATE DATA TO SITE LEVEL, FILTER AND SAVE MODEL READY DATASETS  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;06.0_finalize_data_aggregate_site_level.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;06.1_finalize_data_filter_site_level.R  
07 EXPLORATORY DATA ANALYSIS  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;07.0_exploratory_response_variable_distribution.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;07.1_exploratory_pairs_site_level.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;07.2_exploratory_spatial_blocking.R  
08 RUN MODELS - CROSS VALIDATION  (HPC)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;08.0_nested_spatial_cv_binary_loocv.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;08.1_nested_spatial_cv_continuous_loocv.R  
09 CROSS-VALIDATION FIGURES AND MODEL COMPARISON  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;09.0_cv_model_accuracy_binary.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;09.1_cv_model_accuracy_continuous.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;09.2_cv_tuning_evaluation_binary.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;09.3_cv_tuning_evaluation_continuous.R  
10 CROSS-VALIDATION FINAL MODEL METRICS AND FIGURES  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10.0_cv_final_model_accuracy_plot.R  
11 RUN MODELS - FINAL MONTE CARLO (HPC)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;11.0_mc_final_models_binary.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;11.1_mc_final_models_continuous.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;11.2_mc_final_models_binary_thresholds.R  
12 FORMAT MODELS FOR GEE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;12.0_mc_format_model_for_gee_mc.R  
13 ENSURE R AND GEE PREDICTIONS MATCH  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;13.0_mc_compare_r_gee_predictions.R  
14 VARIABLE IMPORTANCE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;14.0_mc_vip.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;14.1_mc_pdp.R  
16 BIOMASS SUMMARIES  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;16.0_biomass_summaries.R  
17 FIGURES  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;17.0_figures_zone_map.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;17.1_figures_vip.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;17.2_figures_biomass_density_cavm.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;17.3_figures_hillshade_3d.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;17.4_figures_topography_transect.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;17.5_figures_climate_analysis.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;17.6_figures_compare_pixels.R  
18 FINAL MAPS  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;18.0_tile_index_map.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;18.1_compress_final_maps.R  

### Google Earth Engine code

All GEE scripts are Python (earthengine-api).

**utils/** (shared function package)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ccdc.py  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;field_data.py  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;landsat.py  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;misc.py  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;palettes.py  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;params.py  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;pipeline.py  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;refl.py  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sentinel.py  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;temporal_segmentation.py  

**full_workflow/** (pan-Arctic biomass mapping pipeline)  
000 GUIDES AND TEMPLATES  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;000_script_guide  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;000_script_template  
01 PROCESS CAL/VAL DATA  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;01_0_data_unveg_random_pts  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;01_1_data_site_representativeness  
02 PREPARE DATA FOR CCDC MODELING  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;02_0_prep_finalize_arctic_roi  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;02_1_prep_create_water_ice_masks  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;02_2_prep_group_tiles  
03 RUN CCDC  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;03_0_ccdc_fit  
04 PROCESS SEASONAL REFLECTANCE DATA  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;04_0_proc_seasonal_doys  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;04_1_proc_seasonal_doys_median  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;04_2_proc_seasonal_reflectance  
05 TOPOGRAPHIC CORRECTION  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;05_0_topocorr_seasonal_doy_min_max  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;05_1_topocorr_calc_copernicus_dem_slope  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;05_2_topocorr_calc_illumination_condition  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;05_3_topocorr_apply  
06 PREDICTORS  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;06_0_predictor_ecoregion  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;06_1_predictor_wte  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;06_2_predictor_aux_refl  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;06_3_predictor_permafrost  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;06_4_predictor_terraclimate  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;06_5_predictor_tree_cover  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;06_6_predictor_zone  
07 PREDICTOR DATA EXTRACTION  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;07_0_extract_calval  
08 RUN MODELS  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;08_0_run_model_compare_gee_r  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;08_1_run_model_binary_uncertainty  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;08_2_run_model_continuous_uncertainty  
09 EXPORT DATA  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;09_0_export_maps_to_asset  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;09_1_export_biomass_mask  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;09_2_export_mc_iterations_to_gcs  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;09_3_export_final_maps_to_drive  
10 ANALYSES  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10_0_analysis_summarize_biomass  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10_1_analysis_climate_sum_area  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10_2_analysis_climate_mean_biomass  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10_3_analysis_climate_mean_cavm  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10_4_analysis_compare_pixels  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10_5_analysis_masked_area  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10_6_analysis_zone_area  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10_7_analysis_soil_carbon  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10_8_analysis_topography_gradient  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;10_9_analysis_hillshade_data  

**annual/** (annual biomass time series pipeline)  
00 PRE-PROCESSING (run once per model version)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;00_0_modfiles_ecoregion_lookup  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;00_1_modfiles_shorten_strings  
01 SEASONAL REFLECTANCE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;01_0_run_annual_reflectance  
01 MODEL INFERENCE  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;01_1_run_annual_binary_continuous  
01 FINAL BIOMASS MAPS  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;01_2_run_annual_final_biomass  

