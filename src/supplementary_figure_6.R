###################################
###   SUPPLEMENTARY_FIGURE_6   ###
#################################

# Generation of Supplementary Figure 6

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

# check sample categories
table(GDSC_age$status)
table(GDSC_age$tissue_or_cancer_chosen)
# plot sample categories
plot_sample_categories_pan_drug(GDSC_age, 10, 
                                output_dir ="Supplementary Figures/Supplementary_Figure_6.pdf")
