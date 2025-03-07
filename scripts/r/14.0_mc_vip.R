################################################################################
################################################################################

# DESCRIPTION:

# This script summarizes predictor importance

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

in_dir = 'output/11_mc_final_models/'
out_dir = 'output/14_mc_vip/'

# 1.3 Read in data -------------------------------------------------------------

# Permutation 

vip_permutation_files_total_binary = list.files(paste0(in_dir, 'binary/vip/'), paste0('vip_permutation_total_binary'), full.names = TRUE)
vip_permutation_files_woody_binary = list.files(paste0(in_dir, 'binary/vip/'), paste0('vip_permutation_woody_binary'), full.names = TRUE)
vip_permutation_files_total_continuous = list.files(paste0(in_dir, 'continuous/vip/'), paste0('vip_permutation_total_continuous'), full.names = TRUE)
vip_permutation_files_woody_continuous = list.files(paste0(in_dir, 'continuous/vip/'), paste0('vip_permutation_woody_continuous'), full.names = TRUE)

vip_permutation_total_binary = data.frame()
for(file in vip_permutation_files_total_binary){
  mc_num = as.numeric(tail(strsplit(strsplit(file, '[.]')[[1]][1], '_')[[1]], 1))
  df = read.csv(file)
  df$mc = mc_num
  vip_permutation_total_binary = bind_rows(vip_permutation_total_binary, df)
}

vip_permutation_woody_binary = data.frame()
for(file in vip_permutation_files_woody_binary){
  mc_num = as.numeric(tail(strsplit(strsplit(file, '[.]')[[1]][1], '_')[[1]], 1))
  df = read.csv(file)
  df$mc = mc_num
  vip_permutation_woody_binary = bind_rows(vip_permutation_woody_binary, df)
}

vip_permutation_total_continuous = data.frame()
for(file in vip_permutation_files_total_continuous){
  mc_num = as.numeric(tail(strsplit(strsplit(file, '[.]')[[1]][1], '_')[[1]], 1))
  df = read.csv(file)
  df$mc = mc_num
  vip_permutation_total_continuous = bind_rows(vip_permutation_total_continuous, df)
}

vip_permutation_woody_continuous = data.frame()
for(file in vip_permutation_files_woody_continuous){
  mc_num = as.numeric(tail(strsplit(strsplit(file, '[.]')[[1]][1], '_')[[1]], 1))
  df = read.csv(file)
  df$mc = mc_num
  vip_permutation_woody_continuous = bind_rows(vip_permutation_woody_continuous, df)
}

# SHAP 

vip_shap_files_total_binary = list.files(paste0(in_dir, 'binary/vip/'), paste0('vip_shap_total_binary'), full.names = TRUE)
vip_shap_files_woody_binary = list.files(paste0(in_dir, 'binary/vip/'), paste0('vip_shap_woody_binary'), full.names = TRUE)
vip_shap_files_total_continuous = list.files(paste0(in_dir, 'continuous/vip/'), paste0('vip_shap_total_continuous'), full.names = TRUE)
vip_shap_files_woody_continuous = list.files(paste0(in_dir, 'continuous/vip/'), paste0('vip_shap_woody_continuous'), full.names = TRUE)

vip_shap_total_binary = data.frame()
for(file in vip_shap_files_total_binary){
  mc_num = as.numeric(tail(strsplit(strsplit(file, '[.]')[[1]][1], '_')[[1]], 1))
  df = read.csv(file)
  df$mc = mc_num
  vip_shap_total_binary = bind_rows(vip_shap_total_binary, df)
}

vip_shap_woody_binary = data.frame()
for(file in vip_shap_files_woody_binary){
  mc_num = as.numeric(tail(strsplit(strsplit(file, '[.]')[[1]][1], '_')[[1]], 1))
  df = read.csv(file)
  df$mc = mc_num
  vip_shap_woody_binary = bind_rows(vip_shap_woody_binary, df)
}

vip_shap_total_continuous = data.frame()
for(file in vip_shap_files_total_continuous){
  mc_num = as.numeric(tail(strsplit(strsplit(file, '[.]')[[1]][1], '_')[[1]], 1))
  df = read.csv(file)
  df$mc = mc_num
  vip_shap_total_continuous = bind_rows(vip_shap_total_continuous, df)
}

vip_shap_woody_continuous = data.frame()
for(file in vip_shap_files_woody_continuous){
  mc_num = as.numeric(tail(strsplit(strsplit(file, '[.]')[[1]][1], '_')[[1]], 1))
  df = read.csv(file)
  df$mc = mc_num
  vip_shap_woody_continuous = bind_rows(vip_shap_woody_continuous, df)
}

################################################################################
################################################################################

# ==============================================================================
# 2. TIDY ======================================================================
# ==============================================================================

# 2.1 Combine ------------------------------------------------------------------

# Permutation
vip_permutation_total_binary = vip_permutation_total_binary %>% mutate(ds_type = 'Plant', response_type = 'Classification', importance_type = 'Permutation')
vip_permutation_woody_binary = vip_permutation_woody_binary %>% mutate(ds_type = 'Woody Plant', response_type = 'Classification', importance_type = 'Permutation')
vip_permutation_total_continuous = vip_permutation_total_continuous %>% mutate(ds_type = 'Plant', response_type = 'Regression', importance_type = 'Permutation')
vip_permutation_woody_continuous = vip_permutation_woody_continuous %>% mutate(ds_type = 'Woody Plant', response_type = 'Regression', importance_type = 'Permutation')

# SHAP
vip_shap_total_binary = vip_shap_total_binary %>% mutate(ds_type = 'Plant', response_type = 'Classification', importance_type = 'SHAP')
vip_shap_woody_binary = vip_shap_woody_binary %>% mutate(ds_type = 'Woody Plant', response_type = 'Classification', importance_type = 'SHAP')
vip_shap_total_continuous = vip_shap_total_continuous %>% mutate(ds_type = 'Plant', response_type = 'Regression', importance_type = 'SHAP')
vip_shap_woody_continuous = vip_shap_woody_continuous %>% mutate(ds_type = 'Woody Plant', response_type = 'Regression', importance_type = 'SHAP')

vip_all = rbind(vip_permutation_total_binary, vip_permutation_woody_binary, vip_permutation_total_continuous, vip_permutation_woody_continuous,
                vip_shap_total_binary, vip_shap_woody_binary, vip_shap_total_continuous, vip_shap_woody_continuous)

# 2.2 Add count ----------------------------------------------------------------

# For each predictor, denotes how many times (out of 100 Monte Carlo iterations)
# it was selected

vip_all = vip_all %>%
  group_by(ds_type, response_type, importance_type) %>%
  add_count(Variable)

# 2.3 Find and populate missing cases ------------------------------------------

# Complete missing cases
vip_all = data.frame(vip_all) %>%
  complete(Variable, mc, ds_type, response_type, importance_type)

# Where NA, change to 0
vip_all$Importance[is.na(vip_all$Importance)] = 0
vip_all$n[is.na(vip_all$n)] = 0

# To get accurate counts, group and take max of n
vip_all = vip_all %>%
  group_by(ds_type, response_type, importance_type, Variable) %>%
  mutate(n = max(n))

################################################################################
################################################################################

# ==============================================================================
# 3. SUMMARIZE =================================================================
# ==============================================================================

# 3.1 Average predictor importance across Monte Carlo iterations ---------------

# Average feature importance across MC iterations
# Divide by total number of MC iterations, assume if not chosen for iteration importance = 0

vip_avg = vip_all %>%
  group_by(ds_type, response_type, importance_type, Variable) %>%
  summarize(Importance_avg = if_else(importance_type == 'SHAP',
                                     mean(abs(Importance)), # For SHAP, take absolute value
                                     mean(Importance)), # For permutation, do NOT take absolute value
            n = first(n),
            Importance_proportion = n/100) %>%
  distinct()

# 3.2 Calculate categorical predictor importance -------------------------------

# Summing across categories as recommended here:
# https://github.com/shap/shap/issues/397

vip_zone = vip_avg %>%
  group_by(ds_type, response_type, importance_type) %>%
  filter(grepl('zone_3', Variable)) %>%
  summarize(Importance_avg = sum(Importance_avg),
            n = mean(n),
            Importance_proportion = mean(Importance_proportion)) %>%
  mutate(Variable = 'bioclimate zone')

vip_ecosystem = vip_avg %>%
  group_by(ds_type, response_type, importance_type) %>%
  filter(grepl('world_terrestrial_ecosystems', Variable)) %>%
  summarize(Importance_avg = sum(Importance_avg),
            n = mean(n),
            Importance_proportion = mean(Importance_proportion)) %>%
  mutate(Variable = 'ecosystem')

vip_final = vip_avg %>%
  filter(!grepl('zone_3|world_terrestrial_ecosystems', Variable)) %>%
  bind_rows(vip_zone, vip_ecosystem)

# 3.3 Get top predictors -------------------------------------------------------

vip_topN_value = vip_final %>%
  group_by(ds_type, response_type, importance_type) %>%
  slice_max(Importance_avg, n = top_n)

vip_topN_prop = vip_final %>%
  group_by(ds_type, response_type, importance_type) %>%
  slice_max(Importance_proportion, n = top_n)
vip_topN_prop = vip_topN_prop %>% filter(importance_type == 'Permutation')
vip_topN_prop$Importance_avg = vip_topN_prop$Importance_proportion
vip_topN_prop$importance_type = 'Proportion'

vip_topN = rbind(vip_topN_value, vip_topN_prop)

# 3.4 Scale predictor importance -------------------------------------------------------

# For regression model permutation importance units are mean squared error
# Need to take square root to get RMSE and into units of response

vip_topN[vip_topN$importance_type == 'Permutation' & vip_topN$response_type == 'Regression', ]$Importance_avg = sqrt(vip_topN[vip_topN$importance_type == 'Permutation' & vip_topN$response_type == 'Regression', ]$Importance_avg)

################################################################################
################################################################################

# ==============================================================================
# 4. SAVE ======================================================================
# ==============================================================================

write.csv(vip_topN, paste0(out_dir, 'vip_top', top_n, '.csv'), row.names = FALSE)

################################################################################
################################################################################

# ==============================================================================
# 5. PLOT ======================================================================
# ==============================================================================

# 5.1 Tidy for plotting --------------------------------------------------------

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

# 5.2 Plot --------------------------------------------------------

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