###################################
###   SUPPLEMENTARY_FIGURE_19   ###
#################################

# Generation of Supplementary Figure 19

# ============================

# set working directory
setwd("/group/iorio/Alessandro.D/epiclock")

library(ggplot2)

source("src/global_params.R")
source("src/plot_generation_helpers.R")

# import data
df_in <- read.csv(paste0(results_path, "cell_lines/predictions/cls_age_at_sampling_vs_age_prediction.csv"))

# plot
plot_age_scatter_with_compact_legend(
  df = df_in,
  output_path = "Supplementary Figures/Supplementary_Figure_19.pdf",
  legend_position = "right",
  legend_ncol = 2,
  width = 12,
  height = 9,
  stats_text_size = 4
)
