################################################################################
################################################################################

# DESCRIPTION:

# This script creates figures showing predictor importance

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

# NOTES:

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpmisc)
library(viridis)
library(probably)
library(grid)
library(stringr)
library(scales)
library(ggh4x)
library(tidytext)

# 1.2 Parameters ---------------------------------------------------------------

top_n = 10
plt_vip_type = 'Permutation' # Choose 'Permutation', 'SHAP', or 'Proportion'
plt_base_size = 50
plt_eq_size = plt_base_size/3

in_dir = 'output/14_mc_vip/'
out_dir = 'output/17_figures'

# 1.3 Read in data -------------------------------------------------------------

vip_topN = read.csv(paste0(in_dir, 'vip_top', top_n, '.csv'), row.names = FALSE)

# ==============================================================================
# 2. PLOT ======================================================================
# ==============================================================================

# 2.1 Tidy for plotting --------------------------------------------------------

# Tidy predictor names
vip_topN$Variable = gsub('spectral_', '', vip_topN$Variable)
vip_topN$Variable = gsub('_', ' ', vip_topN$Variable)
vip_topN$Variable = gsub('start', '(start ', vip_topN$Variable)
vip_topN$Variable = gsub('early', '(early ', vip_topN$Variable)
vip_topN$Variable = gsub('peak', '(peak ', vip_topN$Variable)
vip_topN$Variable = gsub('late', '(late ', vip_topN$Variable)
vip_topN$Variable = gsub('end', '(end ', vip_topN$Variable)
vip_topN$Variable = gsub('annual', '(annual ', vip_topN$Variable)
vip_topN$Variable = gsub('change', '(change ', vip_topN$Variable)
vip_topN$Variable = gsub('Snowfree', 'snowfree)', vip_topN$Variable)
vip_topN$Variable = gsub('Summer', 'summer)', vip_topN$Variable)
vip_topN$Variable = gsub('Median', 'median)', vip_topN$Variable)
vip_topN$Variable = gsub('Mean', 'mean)', vip_topN$Variable)
vip_topN$Variable = gsub('Range', 'range)', vip_topN$Variable)
vip_topN$Variable = gsub('SSES', 'SSES)', vip_topN$Variable)
vip_topN$Variable = gsub('ESPS', 'ESPS)', vip_topN$Variable)
vip_topN$Variable = gsub('PSLS', 'PSLS)', vip_topN$Variable)
vip_topN$Variable = gsub('LSES', 'LSES)', vip_topN$Variable)
vip_topN$Variable = gsub('trees', 'tree', vip_topN$Variable)
vip_topN$Variable = gsub(' tpi', ' topographic position index', vip_topN$Variable)
vip_topN$Variable = gsub(' cti', ' compound topographic index', vip_topN$Variable)
vip_topN$Variable = gsub(' dem', ' elevation', vip_topN$Variable)
vip_topN$Variable = gsub('topo ', '', vip_topN$Variable)

# Re-order factor
vip_topN$importance_type = factor(vip_topN$importance_type, levels=c('Proportion', 'Permutation', 'SHAP'))

# 2.2 Plot --------------------------------------------------------

plt = ggplot(data = vip_topN[vip_topN$importance_type == plt_vip_type,], aes(y = reorder_within(Variable, Importance_avg, list(importance_type, response_type, ds_type)), x = Importance_avg))+
  ggh4x::facet_grid2(response_type ~ ds_type, scales = "free", independent = "all")+
  geom_col()+
  geom_text(data = vip_topN[vip_topN$importance_type == plt_vip_type,], aes(label = paste0(round(n, 0), '%')), size = 8, hjust = -0.25)+ 
  theme_minimal(base_size = plt_base_size)+
  scale_x_continuous(labels = label_comma(), 
                     expand = expansion(mult = c(0, 0.2)),
                     limits = c(0, NA))+
  scale_y_reordered()+
  labs(y = '', x = 'Variable Importance')+
  theme(axis.text.x = element_text(size = plt_base_size/2),
        axis.text.y = element_text(size = plt_base_size/2),
        plot.margin = unit(c(0,0,0,0), "cm"))
plt

# 5.3 Save --------------------------------------------------------

ggsave(
  paste0(out_dir, 'vip_', tolower(plt_vip_type), '_top', top_n, '.jpg'),
  plt,
  width = 55,
  height = 40,
  units = 'cm',
  bg = 'white',
  dpi = 600
)