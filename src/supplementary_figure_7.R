###################################
###   SUPPLEMENTARY_FIGURE_7   ###
#################################

# Generation of Supplementary Figure 7

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

# load GDSC x cell lines cancer specific data (filtered after category selection)
GDSC_age_filtered <- readRDS(paste0(results_path,
                                    "cell_lines/cell_lines_drug_sensitivity/cell_lines_per_drugs_cancer_specific.rds"))

# plot sample categories
plot_sample_categories_cancer_specific_drug(GDSC_age_filtered,
                                            output_dir ="Supplementary Figures/Supplementary_Figure_7.pdf" )
