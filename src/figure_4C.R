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

# load GDSC x cell lines cancer specific data (filtered after category selection)
GDSC_age_filtered <- readRDS(paste0(results_path,
                                     "cell_lines/cell_lines_drug_sensitivity/cell_lines_per_drugs_cancer_specific.rds"))

# load anova cumulative cancer-specific correlation results
cumulative_correlation_results<- read_xlsx(paste0(results_path,
                                                  "cell_lines/cell_lines_drug_sensitivity/anova_cancer_specific_cumulative_correlation_results.xlsx"))

# anova results
cancer_specific_anova_results <- anova_drugresponse_cancer_specific(GDSC_age_filtered, "age_prediction")

# extract list of significant drugs per cancer type
significant_drugs_list <- extract_significant_drugs(cancer_specific_anova_results)
print(significant_drugs_list)

# plot volcano plots highlighting significant drug-age correlations in these cancer types
plot_cancer_specific_volcano(cumulative_correlation_results, c("Bladder Carcinoma",
                                                               
                                                               "Head and Neck Carcinoma",
                                                               "Oral Cavity Carcinoma"),
                             output_dir = "Figures/Figure4C1.pdf")

# individual correlation plots for statistically significant drugs
plot_significant_drugs(GDSC_age_filtered, significant_drugs_list,
                       output_dir = "Figures/")

