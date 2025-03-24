####################################################
###   ANOVA drug response                        ###
###   Author: Digilio Alessandro                 ###
###   Date: 05/03/2025                           ###
####################################################

# in this script is performed the correlation between the drug sensitivity of the cell lines 
# ln(IC50) data from GDSC (corrected with ANOVA from the variability of tissue/cancer type and 
# MSI status) and the predicted age from the new model

library(tidyverse)
library(dplyr)
library(data.table)
library(fgsea)
library(writexl)
library(MutExMatSorting)
library(reshape2)

# set working directory
setwd("/group/iorio/Alessandro.D/epiclock")

# load scripts
source("src/global_params.R")
source("src/utils.R")
source("src/plot_generation_helpers.R")
source("src/load_data.R")

# set seed for reproducibility
set.seed(123)

# ============================
# DATA LOADING & PREPROCESSING
# ============================

# merge cell line annotations and drug sensitivity data
# applying the selection algorithm for categories selection
GDSC_age <- generate_merged_drug_dataset(GDSC_combined, cl_samples, threshold = 10)
head(GDSC_age)
print(paste("In the pan-cancer analysis, the number of cell lines after the category selection is", length(unique(GDSC_age$COSMIC_ID))))

# check sample categories
table(GDSC_age$status)
table(GDSC_age$tissue_or_cancer_chosen)
# plot sample categories
plot_sample_categories_pan_drug(GDSC_age, 10)
plot_included_tissue_pan_drug(GDSC_age)

# ============================
# PAN-CANCER ANALYSIS
# ============================

# ANOVA MODEL FOR DRUG RESPONSE (Pan-Cancer)
# ------------------------------------------
# here is performed an ANOVA analysis to evaluate how tissue/cancer type and MSI status influence drug sensitivity (LN_IC50) across cell lines.
# The model accounts for the effects of these factors and extracts residuals.Residuals represent the unexplained variance in drug response after 
# considering tissue/cancer type and MSI status.A correlation analysis is performed between residuals and predicted age to assess whether age 
# influences drug response.False Discovery Rate (FDR) correction is applied to control for multiple comparisons, identifying drugs with significant 
# age-response associations.

# perform anova anova model for drug response
anova_drug_pancancer_results <- anova_analysis_drugresponse_pan(GDSC_age)
# plot correlation between predicted age and residual drug sensitivity
plot_predage_LNIC50_correlation_pan_drug(anova_drug_pancancer_results$correlation_results)
# check normality of residuals from ANOVA model
plot_anova_LNIC50_residuals(anova_drug_pancancer_results$residuals)

# DRUG SET ENRICHMENT ANALYSIS (DSEA) - PAN-CANCER
# ------------------------------------------------
# perform Drug Set Enrichment Analysis (DSEA) using drugs ranked by their correlation 
# between predicted age and residual drug sensitivity (LN_IC50).
# identifies enriched drug targets associated with age-dependent drug responses.
dsea_pan_cancer_results <- run_dsea_pan_cancer(anova_drug_pancancer_results$correlation_results, GDSC_age)

# plot drug target enrichment and associated drugs
ordered_targets <- dsea_pancancer_dotplot(dsea_result = dsea_pan_cancer_results, "Drug Target Enrichment by Predicted Group")
plot_target_drug_distribution(dsea_pan_cancer_results, ordered_targets)

# ============================
# CANCER TYPE-SPECIFIC ANALYSIS
# ============================

# filter dataset to include only cell lines categorized as cancer types with sufficient samples (≥10)
GDSC_age_filtered <- GDSC_age[GDSC_age$cancer_specific == "Included as Cancer Type", ]
GDSC_age_filtered$residuals <- NULL
head(GDSC_age_filtered)
print(paste("In the cancer type-specific analysis, the number of cell lines after the category selection is", length(unique(GDSC_age_filtered$COSMIC_ID))))

# plot sample categories
plot_sample_categories_cancer_specific_drug(GDSC_age_filtered)

# ANOVA MODEL FOR DRUG RESPONSE (Cancer-Specific)
# -----------------------------------------------
# performs the same ANOVA model within each cancer type to analyze drug response, extracts residuals
# and assesses correlation with predicted age.
# the aim is the identification of cancer types where drug response is significantly age-associated.

# perform anova anova model for drug response
cancer_specific_anova_results <- anova_drugresponse_cancer_specific(GDSC_age_filtered)
summarize_significant_drugs(cancer_specific_anova_results)

# extract and plot cumulative correlation results across cancer types
cumulative_correlation_results <- cancer_specific_anova_results$cumulative_correlation_results
plot_predage_LNIC50_cumulative_correlation_cancer_specific_drug(cumulative_correlation_results)

# identify cancer types with statistically significant age-drug correlation
# - Bladder Carcinoma
# - Head and Neck Carcinoma
# - Oral Cavity Carcinoma

# extract list of significant drugs per cancer type
significant_drugs_list <- extract_significant_drugs(cancer_specific_anova_results)
# plot volcano plots highlighting significant drug-age correlations in these cancer types
plot_cancer_specific_volcano(cumulative_correlation_results, c("Bladder Carcinoma",
                                                               "Head and Neck Carcinoma",
                                                               "Oral Cavity Carcinoma"))

# individual correlation plots for statistically significant drugs
plot_significant_drugs(GDSC_age_filtered, significant_drugs_list)

# DRUG SET ENRICHMENT ANALYSIS (DSEA) - CANCER-SPECIFIC
# -----------------------------------------------------
# performs DSEA for each cancer type separately.
# ranks drugs based on correlation between predicted age and drug response (LN_IC50).
# identifies enriched drug targets specific to individual cancer types.
dsea_cancer_specific_results <- run_dsea_cancer_specific(cumulative_correlation_results, GDSC_age_filtered)

# plot cancer-type-specific DSEA results
dsea_cancerspecific_dotplot(dsea_cancer_specific_results, "Cancer-Specific DSEA")







