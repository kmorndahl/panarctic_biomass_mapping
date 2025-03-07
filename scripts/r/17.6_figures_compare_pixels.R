################################################################################
################################################################################

# DESCRIPTION:

# This creates a figure comparing the biomass predictions across:
# - Orndahl 30m
# - Spawn 300m
# - Raynolds 8km

# AUTHOR: Kathleen Orndahl
# DATE: 12-4-2024

################################################################################
################################################################################

# ==============================================================================
# 1. SET UP ====================================================================
# ==============================================================================

# 1.1 Packages -----------------------------------------------------------------

library(ggplot2)
library(yardstick)
library(viridis)
library(ggpmisc)
library(stringr)
library(MethComp)

source('scripts/yardstick_bias.R')

# 1.2 Parameters ---------------------------------------------------------------

ds_types = c('raynolds', 'spawn')
textformula = y ~ x
base_size = 34
eq_size = base_size/3

in_dir = 'output/15_gee_output/comparison/'
out_dir = 'output/17_figures/'

# ==============================================================================
# 2. CREATE PLOTTING FUNCTION ==================================================
# ==============================================================================

make_plots = function(ds_type){

  # 2.1 Read in data -----------------------------------------------------------
  
  if(ds_type == 'raynolds'){
    compare = read.csv(paste0(in_dir, ds_type, '_compare_pixels_avg_gm2_perc100.csv'))
  }else if(ds_type == 'spawn'){
    compare = read.csv(paste0(in_dir, ds_type, '_compare_pixels_avg_gm2_perc10_subsample.csv'))
  }else{stop('Aggregation type not recognized')}
  
  # 2.2 Tidy data --------------------------------------------------------------
  
  names(compare) = gsub(ds_type, 'compare', names(compare))
  
  # 2.3 Get biomass max --------------------------------------------------------
  biomass_max = max(max(compare$biomass_density_gm2_compare), max(compare$biomass_density_gm2))
  
  # 2.4 Calculate accuracy metrics ---------------------------------------------
  fit_rmse = round(yardstick::rmse_vec(compare$biomass_density_gm2, compare$biomass_density_gm2_compare), 2)
  fit_mae = round(yardstick::mae_vec(compare$biomass_density_gm2, compare$biomass_density_gm2_compare), 2)
  fit_rsq = round(yardstick::rsq_vec(compare$biomass_density_gm2, compare$biomass_density_gm2_compare), 2)
  fit_relative_rmse = round(fit_rmse/mean(compare$biomass_density_gm2_compare), 2)
  fit_relative_mae = round(fit_mae/mean(compare$biomass_density_gm2_compare), 2)
  fit_bias = round(bias_vec(compare$biomass_density_gm2, compare$biomass_density_gm2_compare), 2)
  
  # 2.5 Total Least Squares calculation ----------------------------------------
  tls = Deming(y=compare$biomass_density_gm2, x=compare$biomass_density_gm2_compare)
  
  # 2.6 Plot -------------------------------------------------------------------
  
  plt_hex = 
    ggplot(compare, aes(x = biomass_density_gm2_compare,  y = biomass_density_gm2))+
    geom_hex(bins = 15)+
    geom_abline(slope = 1, intercept = 0, lty = 2, size = 1.5)+
    geom_abline(slope = tls['Slope'], intercept = tls['Intercept'], col = 'red', size = 3)+
    ylim(c(-1, biomass_max))+
    xlim(c(-1, biomass_max))+
    scale_fill_viridis(option = 'mako', 
                       direction = -1)+
    labs(x = bquote(Biomass~.(str_to_sentence(ds_type))~(g~m^-2)), 
         y = expression(paste('Biomass Orndahl (g ', m^-2, ')')), 
         col = '', 
         fill = '')+
    coord_fixed(ratio = 1)+
    theme_minimal(base_size = base_size)+
    theme(legend.key.height = unit(1.25, "cm"),
          legend.key.width  = unit(0.75, "cm"))+
    stat_poly_eq(size = eq_size, formula = textformula, geom = 'label', aes(label = paste(after_stat(rr.label), sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1.05, vjust = -0.5, color = 'grey20', fill = 'white', label.size = NA, alpha = 0, label.padding = unit(0.01, "lines"))+
    geom_text(size = eq_size, color = 'grey20', data = compare[1,], aes(label = paste0('Bias = ', fit_bias)), x = Inf, y = -Inf, hjust = 1.05, vjust = 0)
  
  plts = list()
  plts[[ds_type]] = plt_hex
  
  return(plts)

}

# ==============================================================================
# 3. APPLY PLOTTING FUNCTION ===================================================
# ==============================================================================

plts_all = lapply(ds_types, make_plots)
plts_all = unlist(plts_all, recursive = FALSE)

# ==============================================================================
# 4. COMBINE PLOTS =============================================================
# ==============================================================================

# Blank space hack: https://github.com/wilkelab/cowplot/issues/31

final = cowplot::plot_grid(NULL,
                           plts_all[['raynolds']] +
                               theme(plot.margin = unit(c(1, 1, 1, 1), "cm")),
                           plts_all[['spawn']] +
                             theme(plot.margin = unit(c(1, 1, 1, 1), "cm")),
                           rel_heights = c(0.2, 1, 1),
                           nrow = 3,
                           labels = c('', 'a)', 'b)'),
                           label_size = base_size + 10,
                           vjust = c(0, -1 -0.01),
                           hjust = c(0, -11, -10.5))

# ==============================================================================
# 5. SAVE ======================================================================
# ==============================================================================

ggsave(
  paste0(out_dir, '/compare_pixels_tls.jpg'),
  final,
  width = 50,
  height = 40,
  units = 'cm',
  bg = 'white',
  dpi = 600
)