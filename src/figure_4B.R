############################
###   FIGURE_4B         ###
###########################

# Generation of Figure 4B

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
dsea_cancer_specific_results <- read_xlsx(paste0(results_path,
                                            "cell_lines/cell_lines_drug_sensitivity/dsea_per_cancer_type/significant_enrichment_dsea_cancerspecific.xlsx"))

# plot figure 4B (drug target enrichment - cancer specific)
# plot cancer-type-specific DSEA results
dsea_cancerspecific_dotplot(dsea_cancer_specific_results, "Cancer-Specific DSEA",
                            output_dir = "Figures/Figure4B.pdf")
