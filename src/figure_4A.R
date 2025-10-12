############################
###   FIGURE_4A         ###
###########################

# Generation of Figure 4A (4A1 and 4A2) 

# ============================
library(tidyverse)
library(dplyr)
library(data.table)
library(fgsea)
library(writexl)
library(MutExMatSorting)
library(reshape2)
library(readxl)

# set working directory
setwd("/group/iorio/Alessandro.D/epiclock")
# load scripts
source("src/global_params.R")
source("src/utils.R")
source("src/plot_generation_helpers.R")

# load enrichment data

dsea_pan_cancer_results <- read_xlsx(paste0(results_path,
                                            "cell_lines/cell_lines_drug_sensitivity/dsea_pancancer/significant_enrichment_dsea_pancancer.xlsx"))

# plot figure 4A1 (drug target enrichment)
ordered_targets <- dsea_pancancer_dotplot(dsea_result = dsea_pan_cancer_results, "Drug Target Enrichment by Predicted Group",
                                          output_dir = "Figures/Figure4A1.pdf")

# plot figure 4A2 (associated number of drugs)
plot_number_target_drug(
  significant_results_DSEA_df = dsea_pan_cancer_results,
  dsea_pan_cancer_results = dsea_pan_cancer_results,
  ordered_targets = ordered_targets,
  output_dir = "Figures/Figure4A2.pdf"
)

# plot figure 4A2 (associated drugs)
plot_target_drug_distribution(dsea_pan_cancer_results, ordered_targets, 
                              output_dir = "Figures/Figure4A2_alternative.pdf")
