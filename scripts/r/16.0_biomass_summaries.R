################################################################################
################################################################################

# DESCRIPTION:

# This script exports biomass summaries across:
# - bioclimate zones
# - CAVM vegetation community types

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

library(tidyverse)
library(ggpattern)
library(ggforce)
library(see)
library(cowplot)

# 1.2 Parameters ---------------------------------------------------------------

version = 'v20240514'
type = 'compare_spawn_masked' # Choose 'zone', 'cavm_fine', 'cavm_coarse', 'compare_raynolds_masked', 'compare_spawn_masked'
tree_mask = 5 # Tree height threshold (m), choose 5 or 200 (no mask)
woody_percent_method = 'woody_percent' # Choose 'woody_percent', 'woody_percent_vegetated_mean' or 'woody_percent_mean'
zone_palette = c("#88CCEE", "#CC6677", "#888888", "#117733")
cavm_coarse_palette = c("#888888", "#56B4E9","#E69F00", "#009E73")

in_dir = 'output/15_gee_output/summaries/'
out_dir = 'output/16_biomass_summaries/'

# 1.3 Read in data -------------------------------------------------------------

biomass = read.csv(paste0(in_dir, 'biomass_', type, '_summary_mask', tree_mask, '_', version, '.csv'))

################################################################################
################################################################################

# ==============================================================================
# 2. TIDY DATA =================================================================
# ==============================================================================

# 2.1 Initial clean up ---------------------------------------------------------

# Remove unwanted columns
biomass = biomass %>% dplyr::select(-one_of(c('system.index', '.geo', 'woody_percent_sum', 'woody_percent_vegetated_sum')))
biomass = biomass %>% dplyr::select(-contains('area'))

# Identify p50, median predicted biomass for easier splitting later
names(biomass) = gsub('biomass_g', 'biomass_predicted_g', names(biomass))

# 2.2 Calculate non-woody biomass ----------------------------------------------
# Necessary for stacked bar charts

if(!type %in% c('compare_raynolds_masked', 'compare_spawn_masked', 'compare_raynolds_unmasked', 'compare_spawn_unmasked')){

  # Calculate non-woody biomass - lwr
  biomass$nonwoody_biomass_lwr_g_sum = biomass$total_biomass_lwr_g_sum - biomass$woody_biomass_lwr_g_sum
  biomass$nonwoody_biomass_lwr_g_mean = biomass$total_biomass_lwr_g_mean - biomass$woody_biomass_lwr_g_mean

  # Calculate non-woody biomass - median
  biomass$nonwoody_biomass_predicted_g_sum = biomass$total_biomass_predicted_g_sum - biomass$woody_biomass_predicted_g_sum
  biomass$nonwoody_biomass_predicted_g_mean = biomass$total_biomass_predicted_g_mean - biomass$woody_biomass_predicted_g_mean

  # Calculate non-woody biomass - upr
  biomass$nonwoody_biomass_upr_g_sum = biomass$total_biomass_upr_g_sum - biomass$woody_biomass_upr_g_sum
  biomass$nonwoody_biomass_upr_g_mean = biomass$total_biomass_upr_g_mean - biomass$woody_biomass_upr_g_mean

}

# 2.3 Summary specific tidying -------------------------------------------------

if(type == 'zone'){
  
  # Calculate woody percent
  biomass$woody_percent = (biomass$woody_biomass_predicted_g_sum/biomass$total_biomass_predicted_g_sum) * 100
  
  # Pivot longer
  biomass = biomass %>%
    dplyr::select(-dsl) %>% 
    pivot_longer(!c(FIRST_zone, woody_percent_mean, woody_percent_vegetated_mean, woody_percent), names_to = "type", values_to = "biomass_g")
  
  # Remove Sub Arctic
  biomass = biomass[biomass$FIRST_zone != 'Sub Arctic',]
  
  # Separate
  biomass = biomass %>% separate(type, c("biomass_type", NA, "uncertainty_type", NA, 'summary_type'))
  
  # Pivot wider
  biomass = biomass %>% pivot_wider(names_from = uncertainty_type, values_from = biomass_g)
  
}else if(type == 'cavm_fine'){
  
  # Lookup table
  lookup = tibble(veg_code = c(seq(1,5,1), seq(21,24,1), seq(31,34,1), seq(41,43,1), seq(91,93,1), 99),
                  veg_description = c("Cryptogam, herb barren (B1)", "Cryptogam, barren complex (B2)", "Non-carbonate mountain complex (B3)", "Carbonate mountain complex (B4)", "Cryptogam, barren, dwarf-shrub complex (B5)", "Graminoid, forb, cryptogam tundra (G1)", "Graminoid, prostrate dwarf-shrub, forb, moss tundra (G2)", "Non-tussock sedge, dwarf-shrub, moss tundra (G3)", "Tussock-sedge, dwarf-shrub, moss tundra (G4)", "Prostrate dwarf-shrub, herb, lichen tundra (P1)", "Prostrate/hemi-prostrate dwarf-shrub, lichen tundra (P2)", "Erect dwarf-shrub, moss tundra (S1)", "Low-shrub, moss tundra (S2)", "Sedge/grass, moss wetland complex (W1)", "Sedge, moss, dwarf-shrub wetland complex (W2)", "Sedge, moss, low-shrub wetland complex (W3)", "Fresh water", "Salt water", "Glaciers", "Non-Arctic areas"),
                  veg_category = c("Barrens", "Barrens", "Barrens", "Barrens", "Barrens", "Graminoid tundras", "Graminoid tundras", "Graminoid tundras", "Graminoid tundras", "Shrub tundras", "Shrub tundras", "Shrub tundras", "Shrub tundras", "Wetlands", "Wetlands", "Wetlands", "Non-vegetated", "Non-vegetated", "Non-vegetated", "Non-Arctic areas"))

  # Classify vegetation description by category
  biomass = biomass %>% inner_join(.,lookup) %>% dplyr::select(-veg_desc)

  # Calculate woody percent
  biomass$woody_percent = (biomass$woody_biomass_predicted_g_sum/biomass$total_biomass_predicted_g_sum) * 100

  # Pivot longer
  biomass = biomass %>%
    pivot_longer(!c(veg_code, veg_description, woody_percent_mean, woody_percent_vegetated_mean, woody_percent, veg_category), names_to = "type", values_to = "biomass_g")
  
  # If biomass is less than zero, force to zero
  # This can happen for the nonwoody category if more woody biomass is predicted than total biomass
  biomass$biomass_g[biomass$biomass_g<0] = 0
  
  # Separate
  biomass = biomass %>% separate(type, c("biomass_type", NA, "uncertainty_type", NA, 'summary_type'))

  # Pivot wider
  biomass = biomass %>% pivot_wider(names_from = uncertainty_type, values_from = biomass_g)

  # Reorder vegetation categories
  biomass$veg_description = factor(biomass$veg_description, levels = c("Non-Arctic areas", "Glaciers", "Salt water", "Fresh water",
                                                                       "Cryptogam, barren, dwarf-shrub complex (B5)", "Carbonate mountain complex (B4)", "Non-carbonate mountain complex (B3)", "Cryptogam, barren complex (B2)", "Cryptogam, herb barren (B1)", 
                                                                       "Sedge, moss, low-shrub wetland complex (W3)", "Sedge, moss, dwarf-shrub wetland complex (W2)", "Sedge/grass, moss wetland complex (W1)",
                                                                       "Tussock-sedge, dwarf-shrub, moss tundra (G4)", "Non-tussock sedge, dwarf-shrub, moss tundra (G3)", "Graminoid, prostrate dwarf-shrub, forb, moss tundra (G2)", "Graminoid, forb, cryptogam tundra (G1)",
                                                                       "Low-shrub, moss tundra (S2)", "Erect dwarf-shrub, moss tundra (S1)", "Prostrate/hemi-prostrate dwarf-shrub, lichen tundra (P2)", "Prostrate dwarf-shrub, herb, lichen tundra (P1)"))
                                   
  biomass$veg_category = factor(biomass$veg_category, levels = c('Non-vegetated', 'Barrens', 'Wetlands', 'Graminoid tundras', 'Shrub tundras', 'Non-Arctic areas'))

}else if(type == 'cavm_coarse'){
  
  # Lookup table
  lookup = tibble(veg_code = c(0, 2, 3, 4, 9),
                  veg_category = c("Barrens", "Graminoid tundras", "Shrub tundras", "Wetlands", "Non-Arctic areas"))
  
  # Classify vegetation description by category
  biomass = biomass %>% inner_join(.,lookup)

  # Calculate woody percent
  biomass$woody_percent = (biomass$woody_biomass_predicted_g_sum/biomass$total_biomass_predicted_g_sum) * 100
  
  # Pivot longer
  biomass = biomass %>%
    pivot_longer(!c(veg_code, woody_percent_mean, woody_percent_vegetated_mean, woody_percent, veg_category), names_to = "type", values_to = "biomass_g")
  
  # If biomass is less than zero, force to zero
  # This can happen for the nonwoody category if more woody biomass is predicted than total biomass
  biomass$biomass_g[biomass$biomass_g<0] = 0
  
  # Separate
  biomass = biomass %>% separate(type, c("biomass_type", NA, "uncertainty_type", NA, 'summary_type'))
  
  # Pivot wider
  biomass = biomass %>% pivot_wider(names_from = uncertainty_type, values_from = biomass_g)
  
  # Reorder vegetation categories
  biomass$veg_category = factor(biomass$veg_category, levels = c('Barrens', 'Wetlands', 'Graminoid tundras', 'Shrub tundras', 'Non-Arctic areas'))
  
}else if(type == 'compare_raynolds_masked' | type == 'compare_raynolds_unmasked'){
  
  # Pivot longer
  biomass = biomass %>%
    dplyr::select(-dsl) %>% 
    pivot_longer(!c(FIRST_zone), names_to = "type", values_to = "biomass_g")
  
  # Filter
  biomass = biomass[biomass$FIRST_zone != 'Sub Arctic',]
  biomass = biomass[biomass$FIRST_zone != 'Oro Arctic',]
  
  # Tidy names
  biomass$type = gsub('_extent', 'Extent', biomass$type)
  
  # Separate
  biomass = biomass %>%
    separate(type, c("biomass_type", NA, "uncertainty_type", NA, 'map_type', 'summary_type'))

  # Pivot wider
  biomass = biomass %>% pivot_wider(names_from = uncertainty_type, values_from = biomass_g)

}else if(type == 'compare_spawn_masked' | type == 'compare_spawn_unmasked'){
  
  # Pivot longer
  biomass = biomass %>%
    dplyr::select(-dsl) %>% 
    pivot_longer(!c(FIRST_zone), names_to = "type", values_to = "biomass_g")
  
  # Filter
  biomass = biomass[biomass$FIRST_zone != 'Sub Arctic',]
  biomass = biomass[biomass$FIRST_zone != 'Oro Arctic',]
  
  # Tidy names
  biomass$type = gsub('_extent', 'Extent', biomass$type)
  biomass$type = gsub('agb_g', 'total_biomass_predicted_g', biomass$type)
  biomass$type = gsub('bgb_g', 'bgb_biomass_predicted_g', biomass$type)
  
  # Separate
  biomass = biomass %>%
    separate(type, c("biomass_type", NA, "uncertainty_type", NA, 'map_type', 'summary_type'))
  
  # Pivot wider
  biomass = biomass %>% pivot_wider(names_from = uncertainty_type, values_from = biomass_g)
  
}

# 2.4 Final tidying ------------------------------------------------------------

# Tidy
biomass$biomass_type = factor(biomass$biomass_type, levels=c('nonwoody', 'woody', 'total', 'bgb'))
biomass$biomass_type = plyr::revalue(biomass$biomass_type, c("total" = "Actual Total", "nonwoody" = "Plant", "woody" = "Woody Plant", "bgb" = "Belowground"))

# For mean and standard deviation, divide by Landsat pixel area (900m2) to get mass/m2
biomass[biomass$summary_type == 'mean',]$lwr = biomass[biomass$summary_type == 'mean',]$lwr / 900
biomass[biomass$summary_type == 'mean',]$predicted = biomass[biomass$summary_type == 'mean',]$predicted / 900
biomass[biomass$summary_type == 'mean',]$upr = biomass[biomass$summary_type == 'mean',]$upr / 900

# Calculate other units
biomass$lwr_Tg = biomass$lwr / 1e+12
biomass$predicted_Tg = biomass$predicted / 1e+12
biomass$upr_Tg = biomass$upr / 1e+12

################################################################################
################################################################################

# ==============================================================================
# 3. SAVE ======================================================================
# ==============================================================================

write.csv(biomass, paste0(out_dir, 'final_biomass_summary_', type, '.csv'), row.names = FALSE)
