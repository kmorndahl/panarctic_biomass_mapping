################################################################################
################################################################################

# DESCRIPTION:

# This script creates partial dependence plots

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

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
library(DALEXtra)
library(pdp)
library(ranger)
library(cowplot)

# 1.2 Parameters ---------------------------------------------------------------

version = 'v20240514'
vip_type = 'Permutation'
response_type = 'continuous'
ds_type = 'total'
top_n = 3
if(response_type == 'continuous'){
  response_var = 'biomass_density_gm2'
  response_name = expression(paste('Biomass (g ', m^-2, ')'))
}else if(response_type == 'binary'){
  response_var = 'presence_absence'
  response_name = 'Probability of Presence'
}

in_dir = paste0('output/11_mc_final_models/', response_type, '/')
out_dir = paste('output/14_mc_vip/')

# 1.3 Read in data -------------------------------------------------------------

# Data
data_files = list.files(paste0(in_dir, 'final_data/'), paste0('final_data_', ds_type, '_', response_type, '_', version, '.*csv'), full.names = TRUE)
data_list = lapply(data_files, read.csv)

# Models
mod_files = list.files(paste0(in_dir, 'models/'), paste0(ds_type, '_', response_type, '_', version, '.*rds'), full.names = TRUE)
mod_list = lapply(mod_files, readRDS)

# Most important predictors
predictors_df = read.csv(paste0(out_dir, 'vip_top10.csv'))

################################################################################
################################################################################

# ==============================================================================
# 2. TIDY AND AGGREGATE DATA ===================================================
# ==============================================================================

# 2.1 Filter data --------------------------------------------------------------

# Restrict to VIP type
predictors_df = predictors_df[predictors_df$importance_type == vip_type,]

# Restrict predictors to top N to fit plot better
predictors_df = predictors_df %>%
  group_by(ds_type, importance_type) %>%
  slice_max(Importance_avg, n = top_n, with_ties = FALSE)

# Get predictors as list
if(ds_type == 'total'){
  predictors = predictors_df[predictors_df$ds_type == 'Plant' & predictors_df$importance_type == vip_type,]$Variable
}else if(ds_type == 'woody'){
  predictors = predictors_df[predictors_df$ds_type == 'Woody Plant' & predictors_df$importance_type == vip_type,]$Variable
}else{stop('Dataset type not recognized')}

# 2.2 Get PDP data for each model ----------------------------------------------

pdp_data_all = data.frame()
for(i in 1:100){
  
  print(paste0('Monte Carlo iteration number ', i))
  cat('\n')
  
  # Get Monte Carlo iteration number
  file = mod_files[[i]]
  mc_num = as.numeric(tail(strsplit(strsplit(file, '[.]')[[1]][1], '_')[[1]], 1))
  
  # Get data
  data = data_list[[i]]

  # Order binary response if necessary
  if(response_type == 'binary'){
    data$presence_absence = as.factor(data$presence_absence)
    data$presence_absence = factor(data$presence_absence, levels = c('presence', 'absence'))
  }
  
  # Get model
  mod = mod_list[[i]]
 
  # Tidy categorical predictors
  mc_data_predictors = names(data)
  mc_data_predictors = gsub('world_terrestrial_', '', mc_data_predictors)
  mc_data_predictors = gsub('s_X.*$', '', mc_data_predictors)
  mc_data_predictors = gsub('s_other', '', mc_data_predictors)
  mc_data_predictors = gsub('_3_.*$', '', mc_data_predictors)
  mc_data_predictors = gsub('zone', 'bioclimate zone', mc_data_predictors)
  mc_data_predictors = unique(mc_data_predictors)
  predictors_tidy = predictors
  
  # Remove if predictor not in current dataset
  predictors_tidy = predictors_tidy[(predictors_tidy %in% mc_data_predictors)]

  # Add back in categorical predictors
  if('ecosystem' %in% predictors_tidy){
    predictors_tidy = predictors_tidy[predictors_tidy != 'ecosystem']
    predictors_tidy = c(predictors_tidy, names(data)[grep('world_terrestrial_ecosystems', names(data))])
  }
  if('bioclimate zone' %in% predictors_tidy){
    predictors_tidy = predictors_tidy[predictors_tidy != 'bioclimate zone']
    predictors_tidy = c(predictors_tidy, names(data)[grep('zone_3', names(data))])
  }
  
  if(length(predictors_tidy) > 0){
    
    # Create explainer
    explainer = explain_tidymodels(
      mod,
      data = dplyr::select(data, -!!sym(response_var)),
      y = as.integer(data[[response_var]]),
      verbose = FALSE
    )
    
    # Create PDP
    pdp = model_profile(explainer, variables = predictors_tidy, N = NULL)
    
    # Get PDP data
    pdp_data = pdp$agr_profiles
    pdp_data$response_type = response_type
    pdp_data$vip_type = vip_type
    pdp_data$ds_type = ds_type
    pdp_data$response_var = response_var
    pdp_data$mc = mc_num
    
    # Append
    pdp_data_all = bind_rows(pdp_data_all, pdp_data)
    
  }
  
}

# 2.3 Format and tidy categorical predictors -----------------------------------

# Format ecosystem predictor
if('ecosystem' %in% predictors){
  wte = pdp_data_all[grep('world_terrestrial_ecosystems', pdp_data_all$`_vname_`),]
  wte = wte %>% pivot_wider(names_from = `_x_`, values_from = `_yhat_`)
  wte$`_yhat_` = wte$`1`
  wte$`_x_` = wte$`_vname_`
  wte$`_vname_` = 'Ecosystem'
  wte = wte %>% select(-c('1', '0'))
  wte = na.omit(wte)
  wte = wte %>% mutate(across(everything(), as.character))
  wte$`_x_` = gsub('world_terrestrial_ecosystems_X', '', wte$`_x_`)
  wte$`_x_` = as.factor(as.numeric(wte$`_x_`))
}else{wte = data.frame()}

# Format zone predictor
if('bioclimate zone' %in% predictors){
  zone = pdp_data_all[grep('zone_3', pdp_data_all$`_vname_`),]
  zone = zone %>% pivot_wider(names_from = `_x_`, values_from = `_yhat_`)
  zone$`_yhat_` = zone$`1`
  zone$`_x_` = zone$`_vname_`
  zone$`_vname_` = 'Bioclimate Zone'
  zone = zone %>% select(-c('1', '0'))
  zone = na.omit(zone)
  zone = zone %>% mutate(across(everything(), as.character))
  zone$`_x_` = gsub('zone_3_', '', zone$`_x_`)
  zone$`_x_` = gsub('[.]', ' ', zone$`_x_`)
  zone$`_x_` = gsub('Arctic', '', zone$`_x_`)
}else{zone = data.frame()}

# Format tree presence predictor
if('trees_presence' %in% predictors){
  trees_presence = pdp_data_all[grep('trees_presence', pdp_data_all$`_vname_`),]
  trees_presence$`_x_` = as.factor(trees_presence$`_x_`)
  trees_presence$`_vname_` = gsub('trees_presence', 'Tree Presence', trees_presence$`_vname_`)
  trees_presence$`_x_` = gsub(0, 'Absent', trees_presence$`_x_`)
  trees_presence$`_x_` = gsub(1, 'Present', trees_presence$`_x_`)
}else{trees_presence = data.frame()}

# Remove categorical predictors from main dataset
pdp_data_continuous = pdp_data_all[!grepl('world_terrestrial_ecosystems', pdp_data_all$`_vname_`),]
pdp_data_continuous = pdp_data_continuous[!grepl('zone_3', pdp_data_continuous$`_vname_`),]
pdp_data_continuous = pdp_data_continuous[!grepl('trees_presence', pdp_data_continuous$`_vname_`),]

# 2.4 Format and tidy continuous predictors ------------------------------------

# Scale spectral predictors
pdp_data_continuous[grepl('spectral', pdp_data_continuous$`_vname_`),]$`_x_` = pdp_data_continuous[grepl('spectral', pdp_data_continuous$`_vname_`),]$`_x_`/10000

# Tidy predictor names
pdp_data_continuous$`_vname_` = gsub('spectral_', '', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('_', ' ', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('start', '(start ', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('early', '(early ', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('peak', '(peak ', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('late', '(late ', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('end', '(end ', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('annual', '(annual ', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('change', '(change ', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('Snowfree', 'snowfree)', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('Summer', 'summer)', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('Median', 'median)', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('Mean', 'mean)', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('Range', 'range)', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('SSES', 'SSES)', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('ESPS', 'ESPS)', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('PSLS', 'PSLS)', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('LSES', 'LSES)', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('trees', 'Tree', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('cover', 'Cover', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub(' tpi', ' Topographic Position Index', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub(' cti', ' Compound Topographic Index', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub(' dem', ' Elevation', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('topo ', '', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('permafrost index', 'Permafrost Index', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('eastness', 'Eastness', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('northness', 'Northness', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('slope', 'Slope', pdp_data_continuous$`_vname_`)
pdp_data_continuous$`_vname_` = gsub('texture NDVI', 'NDVI Texture', pdp_data_continuous$`_vname_`)

predictors = gsub('spectral_', '', predictors)
predictors = gsub('_', ' ', predictors)
predictors = gsub('start', '(start ', predictors)
predictors = gsub('early', '(early ', predictors)
predictors = gsub('peak', '(peak ', predictors)
predictors = gsub('late', '(late ', predictors)
predictors = gsub('end', '(end ', predictors)
predictors = gsub('annual', '(annual ', predictors)
predictors = gsub('change', '(change ', predictors)
predictors = gsub('Snowfree', 'snowfree)', predictors)
predictors = gsub('Summer', 'summer)', predictors)
predictors = gsub('Median', 'median)', predictors)
predictors = gsub('Mean', 'mean)', predictors)
predictors = gsub('Range', 'range)', predictors)
predictors = gsub('SSES', 'SSES)', predictors)
predictors = gsub('ESPS', 'ESPS)', predictors)
predictors = gsub('PSLS', 'PSLS)', predictors)
predictors = gsub('LSES', 'LSES)', predictors)
predictors = gsub('trees', 'Tree', predictors)
predictors = gsub('cover', 'Cover', predictors)
predictors = gsub('presence', 'Presence', predictors)
predictors = gsub(' tpi', ' Topographic Position Index', predictors)
predictors = gsub(' cti', ' Compound Topographic Index', predictors)
predictors = gsub(' dem', ' Elevation', predictors)
predictors = gsub('topo ', '', predictors)
predictors = gsub('permafrost index', 'Permafrost Index', predictors)
predictors = gsub('eastness', 'Eastness', predictors)
predictors = gsub('northness', 'Northness', predictors)
predictors = gsub('slope', 'Slope', predictors)
predictors = gsub('texture NDVI', 'NDVI Texture', predictors)
predictors = gsub('bioclimate zone', 'Bioclimate Zone', predictors)
predictors = gsub('ecosystem', 'Ecosystem', predictors)

# 2.5 Format binary response data ----------------------------------------------
# Need to reverse presence probabilities

if(response_type == 'binary'){
  pdp_data_continuous$`_yhat_` = 1 - as.numeric(pdp_data_continuous$`_yhat_`)
  wte$`_yhat_` = 1 - as.numeric(wte$`_yhat_`)
  zone$`_yhat_` = 1 - as.numeric(zone$`_yhat_`)
  trees_presence$`_yhat_` = 1 - as.numeric(trees_presence$`_yhat_`)
  
}

# 2.6 Combine ------------------------------------------------------------------

pdp_data_final = rbind(pdp_data_continuous, wte, zone, trees_presence)

################################################################################
################################################################################

# ==============================================================================
# 3. PLOT ======================================================================
# ==============================================================================

all_plts = list()
for(predictor in predictors){
  
  df = pdp_data_final[pdp_data_final$`_vname` == predictor,]
  
  if(predictor %in% c('Bioclimate Zone', 'Ecosystem', 'Tree Presence')){
    
    if(predictor == 'Ecosystem'){
      
      plt = ggplot(df, aes(x = `_x_`, y = as.numeric(`_yhat_`)))+
        geom_boxplot(outliers = FALSE)+
        geom_jitter(aes(col = as.factor(mc)), alpha = 0.2, size = 1, width = 0.15, height = 0.15)+
        theme_minimal(base_size = 34)+
        theme(legend.position="none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 30))+
        labs(x = '', y = response_name, title = predictor)

    }else{
    
      plt = ggplot(df, aes(x = `_x_`, y = as.numeric(`_yhat_`)))+
        geom_boxplot(outliers = FALSE)+
        geom_jitter(aes(col = as.factor(mc)), alpha = 0.2, size = 4, width = 0.15, height = 0.15)+
        theme_minimal(base_size = 34)+
        theme(legend.position="none",
              plot.title = element_text(size = 30))+
        labs(x = '', y = response_name, title = predictor)

    }
    
  }else{
    
    plt = ggplot(df, aes(x = as.numeric(`_x_`), y = as.numeric(`_yhat_`)))+
      geom_line(aes(col = as.factor(mc)), alpha = 0.2, lwd = 0.75)+
      geom_smooth(col = 'black', lwd = 2)+
      theme_minimal(base_size = 34)+
      theme(legend.position="none",
            plot.title = element_text(size = 30))+
      labs(x = '', y = response_name, title = predictor)

  }
  
  all_plts[[predictor]] = plt
  
}

final_plt = plot_grid(plotlist = all_plts, nrow = 1)

################################################################################
################################################################################

# ==============================================================================
# 4. SAVE ======================================================================
# ==============================================================================

ggsave(
  paste0(out_dir, 'pdp_', response_type, '_', ds_type, '_', tolower(vip_type), '.jpg'),
  final_plt,
  width = 60,
  height = 16,
  units = 'cm',
  bg = 'white',
  dpi = 600
)

