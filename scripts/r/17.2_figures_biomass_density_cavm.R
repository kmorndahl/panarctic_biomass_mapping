################################################################################
################################################################################

# DESCRIPTION:

# This script creates figure showing biomass across CAVM vegetation community types

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

woody_percent_method = 'woody_percent' # Choose 'woody_percent', 'woody_percent_vegetated_mean' or 'woody_percent_mean'
zone_palette = c("#88CCEE", "#CC6677", "#888888", "#117733")
cavm_coarse_palette = c("#888888", "#56B4E9","#E69F00", "#009E73")

in_dir = 'output/16_biomass_summaries/'
out_dir = 'output/17_figures/'

# 1.3 Read in data -------------------------------------------------------------

biomass = read.csv(paste0(in_dir, 'final_biomass_summary_cavm_fine.csv'))

################################################################################
################################################################################

# ==============================================================================
# 2. PLOT ======================================================================
# ==============================================================================

# 2.1 Tidy -------------------------------------------------------------------

# Filter
biomass = biomass[biomass$veg_category != 'Non-Arctic areas' & biomass$veg_category != 'Non-vegetated',]

# Reorder vegetation categories
biomass$veg_description = factor(biomass$veg_description, levels = c("Non-Arctic areas", "Glaciers", "Salt water", "Fresh water",
                                                                     "Cryptogam, barren, dwarf-shrub complex (B5)", "Carbonate mountain complex (B4)", "Non-carbonate mountain complex (B3)", "Cryptogam, barren complex (B2)", "Cryptogam, herb barren (B1)", 
                                                                     "Sedge, moss, low-shrub wetland complex (W3)", "Sedge, moss, dwarf-shrub wetland complex (W2)", "Sedge/grass, moss wetland complex (W1)",
                                                                     "Tussock-sedge, dwarf-shrub, moss tundra (G4)", "Non-tussock sedge, dwarf-shrub, moss tundra (G3)", "Graminoid, prostrate dwarf-shrub, forb, moss tundra (G2)", "Graminoid, forb, cryptogam tundra (G1)",
                                                                     "Low-shrub, moss tundra (S2)", "Erect dwarf-shrub, moss tundra (S1)", "Prostrate/hemi-prostrate dwarf-shrub, lichen tundra (P2)", "Prostrate dwarf-shrub, herb, lichen tundra (P1)"))
biomass$veg_category = factor(biomass$veg_category, levels = c('Non-vegetated', 'Barrens', 'Wetlands', 'Graminoid tundras', 'Shrub tundras', 'Non-Arctic areas'))

# 2.2 Plot -------------------------------------------------------------------

plt = ggplot(data = biomass[biomass$summary_type == 'mean',], aes(x = veg_description, y = predicted))+
  facet_col(~veg_category, scales = "free_y", space = "free", strip.position = "top")+ # Different pattern spacing
  geom_col_pattern(data = biomass[biomass$summary_type == 'mean' & biomass$biomass_type != 'Actual Total',],
                   aes(pattern = biomass_type, fill = veg_category), 
                   position = position_stack(),
                   colour = "black",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = c(rep(0.04, 10), rep(0.06, 6), rep(0.05, 16)),
                   pattern_key_scale_factor = 0.25,
                   lwd = 0.75)+
  geom_errorbar(data = biomass[biomass$summary_type == 'mean' & biomass$biomass_type == 'Actual Total',],
                aes(ymin=lwr, ymax=upr), 
                width=.2,
                position=position_dodge(.9))+
  scale_pattern_manual(values = c("none", "stripe"),
                       guide = guide_legend(override.aes = list(fill = "white")))+
  scale_fill_manual(values = cavm_coarse_palette)+
  geom_text(aes(label = ifelse(biomass_type == 'Actual Total',  paste0(round(!!sym(woody_percent_method), 0), '%'), ""),
                y = 1750),
            size = 10,
            show.legend = FALSE)+
  labs(y = expression(paste('Biomass (g  ', m^-2, ')')), x = '', col = '', fill = '', pattern = '')+
  guides(color = "none", 
         fill = "none", 
         label = "none",
         pattern = guide_legend(override.aes = list(color = 'black', fill = NA, lwd = 2)))+
  theme_minimal(base_size = 24)+
  theme(legend.position="top", 
        legend.key.width = unit(1, "cm"), 
        legend.key.height = unit(1, "cm"),
        legend.text = element_text(size = 30),
        strip.text = element_text(size = 34, hjust = 0),
        axis.title.x = element_text(size = 34),
        axis.text.x = element_text(size = 28))+
  coord_flip()
  
  ggsave(
    paste0(out_dir, 'biomass_density_cavm.jpg'),
    plt,
    width = 40,
    height = 40,
    units = 'cm',
    bg = 'white',
    dpi = 600
  )



