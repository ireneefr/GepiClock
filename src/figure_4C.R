############################
###   FIGURE_4C         ###
###########################

# Generation of Figure 4C_1 and 4C_specific_drugs

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

# load anova cumulative cancer-specific correlation results
cumulative_correlation_results<- read_xlsx(paste0(results_path,
                                            "cell_lines/cell_lines_drug_sensitivity/anova_cancer_specific_cumulative_correlation_results.xlsx"))

# load significant results
cancer_specific_anova_results <- read_xlsx(paste0(results_path,
                                                  "cell_lines/cell_lines_drug_sensitivity/significant_drugs_cancer_types/significant_drugs_anova_results.xlsx"))

# extract list of significant drugs per cancer type
significant_drugs_list <- extract_significant_drugs(cancer_specific_anova_results)

# plot volcano plots highlighting significant drug-age correlations in these cancer types
plot_cancer_specific_volcano(cumulative_correlation_results, c("Bladder Carcinoma",
                                                               "Head and Neck Carcinoma",
                                                               "Oral Cavity Carcinoma"),
                             output_dir = "Figures/Figure4C1.pdf")

# individual correlation plots for statistically significant drugs
plot_significant_drugs(GDSC_age_filtered, significant_drugs_list,
                       output_dir = "Figures/")
