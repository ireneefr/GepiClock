###################################
###   SUPPLEMENTARY_FIGURE_9   ###
#################################

# Generation of Supplementary Figure 9

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

# load GDSC x cell lines pan cancer data (filtered after category selection)
GDSC_age <- readRDS(paste0(results_path,
                           "cell_lines/cell_lines_drug_sensitivity/cell_lines_per_drugs_pancancer.rds"))
# perform anova anova model for drug response
anova_drug_pancancer_results <- anova_analysis_drugresponse_pan(GDSC_age,"age_prediction")

# check normality of residuals from ANOVA model
plot_anova_LNIC50_residuals(anova_drug_pancancer_results$residuals,
                            output_dir = paste0(figures_path, "anova_LNIC50_residuals_plot.pdf"))

