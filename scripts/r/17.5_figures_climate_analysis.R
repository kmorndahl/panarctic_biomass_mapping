################################################################################
################################################################################

# DESCRIPTION:

# Create figure of climate analysis:
# - Biomass by thawing degree day bin
# - Biomass vs. thawing degree day across CAVM vegetation community types
# - Distribution of thawing degree day bins

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
library(mgcv)
library(pathviewr)

# 1.2 Parameters ---------------------------------------------------------------

gdd_base = 0
correlation_type = 'pearson'
cavm_coarse_palette = c("#888888", "#56B4E9","#E69F00", "#009E73")
plant_woody_plant = c("#E69F00", "#009E73")

if(gdd_base == 0){
  cutoff_filter = 1250
  cor_x = 300
  cor_x_cavm = 300
  x_axis_label = 'Thawing '
}else{
  cutoff_filter = 500
  cor_x = 100
  cor_x_cavm = 100
  x_axis_label = 'Growing '
}

in_dir = 'output/15_gee_output/'
out_dir = 'output/17_figures/'

# 1.3 Read in data -------------------------------------------------------------

biomass_gdd_cavm = read.csv(paste0(in_dir, 'climate/gdd', gdd_base, '/arctic_biomass_gdd', gdd_base, '_mean_1km_cavm_fine.csv'))

sum_dir = paste0(in_dir, 'climate/gdd', gdd_base, '/area_sum/')
sum_files = list.files(sum_dir, paste0('sum.*csv'), full.names = TRUE)
sum_files = sum_files[grep('sample', sum_files, invert = TRUE)]
sum_df = lapply(sum_files, read.csv) %>% bind_rows()

mean_dir = paste0(in_dir, 'climate/gdd', gdd_base, '/biomass_mean')
mean_files = list.files(mean_dir, paste0('mean.*csv'), full.names = TRUE)
mean_files = mean_files[grep('sample', mean_files, invert = TRUE)]
mean_df = lapply(mean_files, read.csv) %>% bind_rows()

################################################################################
################################################################################

# ==============================================================================
# 2. TIDY DATA =================================================================
# ==============================================================================

# 2.1 Combine data -------------------------------------------------------------

biomass_gdd = mean_df
biomass_gdd = biomass_gdd %>% inner_join(sum_df %>% select(c(gdd_bin, area_m2)), by = c('gdd_bin'))

# 2.2 Remove unnecessary columns -----------------------------------------------

biomass_gdd = biomass_gdd %>% select(-c('system.index', '.geo'))
biomass_gdd_cavm = biomass_gdd_cavm %>% select(-c('system.index', '.geo'))

# 2.3 Calculate area percent -----------------------------------------------

total_area = sum(biomass_gdd$area_m2)
biomass_gdd$area_percent = (biomass_gdd$area_m2/total_area)*100

# 2.4 Convert to tidy data format ----------------------------------------------

# Make long
biomass_gdd = biomass_gdd %>% pivot_longer(!any_of(c('gdd_bin', 'area_m2', 'area_percent', 'woody_percent')), names_to = 'ds_type', values_to = 'biomass_gm2')
biomass_gdd_cavm = biomass_gdd_cavm %>% pivot_longer(!any_of(c('gdd', 'gdd_bin', 'area_m2', 'area_percent', 'woody_percent', 'veg_code', 'veg_desc')), names_to = 'ds_type', values_to = 'biomass_gm2')
biomass_gdd$ds_type = gsub('_biomass', '-biomass', biomass_gdd$ds_type)
biomass_gdd = biomass_gdd %>% separate_wider_delim(ds_type, delim = "-", names = c("ds_type", "prediction_type"))
biomass_gdd_cavm$ds_type = gsub('_biomass', '-biomass', biomass_gdd_cavm$ds_type)
biomass_gdd_cavm = biomass_gdd_cavm %>% separate_wider_delim(ds_type, delim = "-", names = c("ds_type", "prediction_type"))

# Make wider
biomass_gdd$prediction_type = gsub('biomass_', '', biomass_gdd$prediction_type)
biomass_gdd$prediction_type = gsub('biomass', 'predicted', biomass_gdd$prediction_type)
biomass_gdd = biomass_gdd %>% pivot_wider(names_from = prediction_type, values_from = biomass_gm2)
biomass_gdd_cavm$prediction_type = gsub('biomass_', '', biomass_gdd_cavm$prediction_type)
biomass_gdd_cavm$prediction_type = gsub('biomass', 'predicted', biomass_gdd_cavm$prediction_type)
biomass_gdd_cavm = biomass_gdd_cavm %>% pivot_wider(names_from = prediction_type, values_from = biomass_gm2)

# 2.5 Format names -------------------------------------------------------------

biomass_gdd$ds_type = gsub('total', 'Plant', biomass_gdd$ds_type)
biomass_gdd$ds_type = gsub('woody', 'Woody Plant', biomass_gdd$ds_type)
biomass_gdd_cavm$ds_type = gsub('total', 'Plant', biomass_gdd_cavm$ds_type)
biomass_gdd_cavm$ds_type = gsub('woody', 'Woody Plant', biomass_gdd_cavm$ds_type)

# 2.6 Calculate relative uncertainty -------------------------------------------

biomass_gdd$relative_uncertainty = (biomass_gdd$upr - biomass_gdd$lwr) / biomass_gdd$predicted

# 2.7 Group and tidy CAVM vegetation community types ---------------------------

# Lookup table
lookup = tibble(veg_code = c(seq(1,5,1), seq(21,24,1), seq(31,34,1), seq(41,43,1), seq(91,93,1), 99),
                veg_code_tidy = c("B1", "B2", "B3", "B4", "B5", "G1", "G2", "G3", "G4", "P1", "P2", "S1", "S2", "W1", "W2", "W3", "Fresh water", "Salt water", "Glaciers", "Non-Arctic areas"),
                veg_category = c("Barrens", "Barrens", "Barrens", "Barrens", "Barrens", "Graminoid tundras", "Graminoid tundras", "Graminoid tundras", "Graminoid tundras", "Shrub tundras", "Shrub tundras", "Shrub tundras", "Shrub tundras", "Wetlands", "Wetlands", "Wetlands", "Non-vegetated", "Non-vegetated", "Non-vegetated", "Non-Arctic areas"))

# Classify vegetation description by category
biomass_gdd_cavm = biomass_gdd_cavm %>% inner_join(.,lookup) %>% dplyr::select(-veg_desc)

# Reorder vegetation categories
biomass_gdd_cavm$veg_code_tidy = factor(biomass_gdd_cavm$veg_code_tidy, levels = c("Non-Arctic areas", "Glaciers", "Salt water", "Fresh water",
                                                                                   "B5", "B4", "B3", "B2", "B1", 
                                                                                   "W3", "W2", "W1",
                                                                                   "G4", "G3", "G2", "G1",
                                                                                   "S2", "S1", "P2", "P1"))

biomass_gdd_cavm$veg_category = factor(biomass_gdd_cavm$veg_category, levels = c('Non-vegetated', 'Barrens', 'Wetlands', 'Graminoid tundras', 'Shrub tundras', 'Non-Arctic areas'))

################################################################################
################################################################################

# ==============================================================================
# 3. PLOT ======================================================================
# ==============================================================================

# 3.1 Determine turning points and GDD cut off ---------------------------------

# Filter data to get monotonically decreasing dataset
turn_df = biomass_gdd[biomass_gdd$ds_type == 'Plant' & biomass_gdd$gdd_bin > cutoff_filter,]
turn_df = turn_df[order(turn_df$gdd_bin),]

# Visualize data
ggplot(turn_df, aes(x = gdd_bin, y = area_m2))+
  geom_point()+
  geom_smooth(method = 'gam')+
  labs(x = 'Growing Degree Days\n(heat sum above 0°C)', y = expression(paste('Total Area in Bin (', m^2, ')')))+
  theme_minimal(base_size = 30)

# Mathematically determine turning point
mod = gam(area_m2 ~ s(gdd_bin, bs = "cs", fx = TRUE, k = 10), data = turn_df)
pred_y = predict(mod)
pred = data.frame(x = turn_df$gdd_bin, y = pred_y)
elbow_idx = find_curve_elbow(pred, plot_curve = TRUE)
elbow = turn_df$gdd_bin[elbow_idx]
elbow_df = data.frame(gdd_bin = turn_df$gdd_bin[elbow_idx], area_m2 = pred_y[elbow_idx])

# Plot elbow
ggplot(turn_df, aes(x = gdd_bin, y = area_m2))+
  geom_point()+
  geom_smooth(method = 'gam')+
  geom_vline(xintercept = elbow, color = 'red')+
  labs(x = 'Growing Degree Days\n(heat sum above 0°C)', y = expression(paste('Total Area in Bin (', m^2, ')')))+
  theme_minimal(base_size = 30)

# Set cutoff
bin_cutoff = elbow

# Calculate area excluded
area_df = biomass_gdd[biomass_gdd$ds_type == 'Plant',]
incl = sum(area_df[area_df$gdd_bin < bin_cutoff,]$area_m2)
excl = sum(area_df[area_df$gdd_bin >= bin_cutoff,]$area_m2)
percent_excl = (excl/(incl+excl))*100
print(paste0('Percentage of data excluded using upper TDD bin threshold: ', round(percent_excl, 2), '%'))

# 3.2 Biomass vs. GDD bin plot -------------------------------------------------

# Calculate correlations
cor_regr = biomass_gdd[biomass_gdd$gdd_bin < bin_cutoff,] %>%
  group_by(ds_type) %>%
  summarise(cor = cor(predicted, gdd_bin, method = correlation_type))
cor_cavm = biomass_gdd_cavm %>%
  filter(ds_type == 'Plant') %>%
  summarise(cor = cor(predicted, gdd, method = correlation_type))

# Plot
plt = ggplot(data = biomass_gdd[biomass_gdd$gdd_bin < bin_cutoff,], aes(x = gdd_bin, y = predicted, col = ds_type))+
  geom_line(lwd = 3)+
  geom_line(aes(x = gdd_bin, y = woody_percent*25), lwd = 3, lty = 'dashed', col = 'grey60')+
  scale_y_continuous(sec.axis = sec_axis(~ ./25, name = "Woody Dominance (%)"))+ 
  geom_ribbon(aes(ymin = lwr, ymax = upr, col = ds_type, fill = ds_type), lwd = 0.00001, alpha = 0.2)+
  scale_color_manual(values = plant_woody_plant)+
  annotate("text", x = cor_x, y = 2200, size = 14, label = paste0('\u03c1 = ', round(cor_regr$cor[cor_regr$ds_type == 'Plant'], 2)), col = "#E69F00")+
  annotate("text", x = cor_x, y = 2000, size = 14, label = paste0('\u03c1 = ', round(cor_regr$cor[cor_regr$ds_type == 'Woody Plant'], 2)), col = "#009E73")+
  theme_minimal(base_size = 60)+
  labs(col = '', x = '', y = expression(paste('Biomass (g  ', m^-2, ')')))+
  guides(fill = "none")+
  theme(legend.position="top", 
        legend.key.width = unit(2, "cm"), 
        legend.key.height = unit(2, "cm"),
        axis.title.y = element_text(),
        axis.title.y.right = element_text(color = 'grey60', size = 50),
        axis.text.y.right = element_text(color = 'grey60'),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-20,-20,-20,-20))+
  guides(color=guide_legend(override.aes=list(fill=NA)))
plt

# Save
ggsave(
  paste0(out_dir, 'biomass_gdd.jpg'),
  plt,
  width = 40,
  height = 40,
  units = 'cm',
  bg = 'white',
  dpi = 600
)

# 3.3 GDD area plot ------------------------------------------------------------

# Plot
plt_area = ggplot(data = biomass_gdd[biomass_gdd$gdd_bin < bin_cutoff & biomass_gdd$ds_type == 'Plant',], 
                  aes(x = gdd_bin, y = area_m2/1000000))+
  geom_col(col = 'grey', fill = 'grey')+
  theme_minimal(base_size = 40)+
  labs(col = '', x = '', y = expression(paste('Area (', km^2, ')')))+
  scale_x_continuous(limits = c(ggplot_build(plt)$layout$panel_params[[1]]$x.range[1], ggplot_build(plt)$layout$panel_params[[1]]$x.range[2]))+
  guides(fill = "none")
plt_area

# Save
ggsave(
  paste0(out_dir, 'gdd_area.jpg'),
  plt_area,
  width = 40,
  height = 40,
  units = 'cm',
  bg = 'white',
  dpi = 600
)

# 3.4 Biomass vs. GDD across vegetation types plot -----------------------------

# Plot
plt_cavm = ggplot(data = biomass_gdd_cavm[!biomass_gdd_cavm$veg_category %in% c('Non-vegetated', 'Non-Arctic areas') & biomass_gdd_cavm$ds_type == 'Plant',], aes(x = gdd, y = predicted, col = veg_category))+
  geom_point(size = 1)+
  geom_errorbar(aes(ymin=lwr, ymax=upr), lwd = 2, lineend = "square", width = 0)+
  geom_label(aes(label = veg_code_tidy, fill = veg_category), color = 'white', size = 6, show.legend = FALSE)+
  theme_minimal(base_size = 40)+
  scale_color_manual(values = cavm_coarse_palette)+
  scale_fill_manual(values = cavm_coarse_palette)+
  annotate("text", x = cor_x_cavm, y = 1500, size = 11, label = paste0('\u03c1 = ', round(cor_cavm$cor, 2)), col = "grey20")+
  theme(legend.text = element_text(size = 20),
        legend.position="top")+
  guides(fill = 'none',
         col = 'none',
         label = 'none',
         text = 'none')+
  labs(col = '', x = '', y = expression(paste('Plant Biomass (g  ', m^-2, ')')), col = '')
plt_cavm

# Save
ggsave(
  paste0(out_dir, 'biomass_gdd_cavm.jpg'),
  plt_cavm,
  width = 40,
  height = 40,
  units = 'cm',
  bg = 'white',
  dpi = 600
)

# 3.5 Final combined plot -----------------------------

# Arrange plots
bottom_row = plot_grid(plt_cavm, plt_area, labels = c('b)', 'c)'), label_size = 60, label_y = 1.175)
plt_all = plot_grid(plt, bottom_row, labels = c('a)', ''), label_size = 60, ncol = 1, rel_heights = c(1.5, 1), label_y = 1)

# Create x axis label
xlabel = ggdraw() + 
  draw_label(paste0(x_axis_label, "Degree Days\n(heat sum above ", gdd_base, "°C)"), x = 0.5, hjust = 0.5, size = 60) +
  theme(plot.margin = margin(0, 0, 0, 0))

# Combine all elements
plt_all_label = plot_grid(plt_all, xlabel, ncol = 1, rel_heights = c(1, 0.1))

# Save
ggsave(
  paste0(out_dir, 'biomass_climate_analysis.jpg'),
  plt_all_label,
  width = 40,
  height = 50,
  units = 'cm',
  bg = 'white',
  dpi = 600
)

