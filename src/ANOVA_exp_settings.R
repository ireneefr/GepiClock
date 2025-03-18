####################################################
###   ANOVA experimental settings                ###
###   Author: Digilio Alessandro                 ###
###   Date: 05/03/2025                           ###
####################################################

# in this script is performed the correlation between some experimental settings data
# (doubling time, mutational burden) of the cell lines (corrected with ANOVA from the 
# variability of tissue/cancer type and  MSI status) and the predicted age 
# from the new model and is compared with the results of the same analysis with donor age 
# at sampling

library(tidyverse)
library(dplyr)
library(data.table)

# set working directory
setwd("/group/iorio/Alessandro.D/EpiClock")

# load scripts
source("src/global_params.R")
source("src/utils.R")
source("src/plot_generation_helpers.R")
source("src/load_data.R")

# set seed for reproducibility
set.seed(123)

# ---------------------------------
# ---------------------------------
# --------- DOUBLING TIME ---------
# ---------------------------------
# ---------------------------------

# load data
# remove samples with NA in doubling_time_hours
cl_samples_dt <- cl_samples_with_age_at_sampling[!is.na(cl_samples_with_age_at_sampling$doubling_time_hours), ]
print(paste("There are", nrow(cl_samples_dt),
            "samples with both 'age at sampling' and 'doubling time' correctly annotated"))

# ============================
# PAN-CANCER ANALYSIS
# ============================

# ANOVA MODEL FOR DOUBLING TIME (Pan-Cancer)
# ------------------------------------------
# here is performed an ANOVA analysis to evaluate how tissue/cancer type and MSI status influence doubling time across cell lines.
# the model accounts for the effects of these factors and extracts residuals. Residuals represent the unexplained variance in doubling time value after 
# considering tissue/cancer type and MSI status. A correlation analysis is performed between residuals and predicted age to assess whether age 
# influences doubling time. The analysis is then performed with "age at sampling" annotation to compare the relationships

# compute correlation between doubling time and predicted age
correlation_dt_pred_age <- anova_exp_settings_analysis(cl_samples_dt, "doubling_time_hours", "age_prediction")

# compute correlation between doubling time and donor age (age at sampling)
correlation_dt_donor_age <- anova_exp_settings_analysis(cl_samples_dt, "doubling_time_hours", "age_at_sampling")

# ============================
# CANCER TYPE-SPECIFIC ANALYSIS
# ============================

# ANOVA MODEL FOR DOUBLING TIME (Cancer-Specific)
# -----------------------------------------------
# performs the same ANOVA model within each cancer type to analyze doubling time, extracts residuals
# and assesses correlation with predicted age.
# the aim is the identification of cancer types where doubling time is significantly age-associated.

# extract samples with the new selection columns
cl_samples_dt <- correlation_dt_pred_age$cl_samples_filtered
cl_samples_dt$residuals <- NULL

# for each cancer type compute correlation between doubling time and predicted age 
correlation_dt_pred_age_cancer_specific <- anova_exp_settings_cancer_specific(cl_samples_dt, "doubling_time_hours", "age_prediction")
# for each cancer type compute correlation between doubling time and donor age (age at sampling)
correlation_dt_donor_age_cancer_specific <- anova_exp_settings_cancer_specific(cl_samples_dt, "doubling_time_hours", "age_at_sampling")

# landscape correlation plot
# extract samples with >= 10 samples per cancer type
cl_samples_dt_filtered <- correlation_dt_pred_age_cancer_specific$cl_samples_filtered
landscape_correlation_plot_dt_age <- plot_landscape_correlation_exp_settings_age(
  cl_samples_dt_filtered, "doubling_time_hours", "age_prediction", "age_at_sampling", cancer_colors
)

# ---------------------------------
# ---------------------------------
# ------- MUTATIONAL BURDEN -------
# ---------------------------------
# ---------------------------------

# load data
# remove samples with NA in mutational_burden
cl_samples_mut <- cl_samples_with_age_at_sampling[!is.na(cl_samples_with_age_at_sampling$mutational_burden), ]
print(paste("There are", nrow(cl_samples_mut),
            "samples with both 'age at sampling' and 'mutational burden' correctly annotated"))

# ============================
# PAN-CANCER ANALYSIS
# ============================

# ANOVA MODEL FOR MUTATIONAL BURDEN (Pan-Cancer)
# ------------------------------------------
# here is performed an ANOVA analysis to evaluate how tissue/cancer type and MSI status influence mutational burden across cell lines.
# the model accounts for the effects of these factors and extracts residuals. Residuals represent the unexplained variance in mutational burden value after 
# considering tissue/cancer type and MSI status. A correlation analysis is performed between residuals and predicted age to assess whether age 
# influences mutational burden. The analysis is then performed with "age at sampling" annotation to compare the relationships

# compute correlation between mutational burden and predicted age
correlation_mut_pred_age <- anova_exp_settings_analysis(cl_samples_mut, "mutational_burden", "age_prediction")

# compute correlation between mutational burden and donor age (age at sampling)
correlation_mut_donor_age <- anova_exp_settings_analysis(cl_samples_mut, "mutational_burden", "age_at_sampling")

# ============================
# CANCER TYPE-SPECIFIC ANALYSIS
# ============================

# ANOVA MODEL FOR MUTATIONAL BURDEN (Cancer-Specific)
# -----------------------------------------------
# performs the same ANOVA model within each cancer type to analyze mutational_burden, extracts residuals
# and assesses correlation with predicted age.
# the aim is the identification of cancer types where mutational_burden is significantly age-associated.

# extract samples with the new selection columns
cl_samples_mut <- correlation_mut_pred_age$cl_samples_filtered
cl_samples_mut$residuals <- NULL

# for each cancer type compute correlation between mutational_burden and predicted age 
correlation_mut_pred_age_cancer_specific <- anova_exp_settings_cancer_specific(cl_samples_mut, "mutational_burden", "age_prediction")
# for each cancer type compute correlation between mutational_burden and donor age (age at sampling)
correlation_mut_donor_age_cancer_specific <- anova_exp_settings_cancer_specific(cl_samples_mut, "mutational_burden", "age_at_sampling")

# landscape correlation plot
# extract samples with >= 10 samples per cancer type
cl_samples_mut_filtered <- correlation_mut_pred_age_cancer_specific$cl_samples_filtered
landscape_correlation_plot_mut_age <- plot_landscape_correlation_exp_settings_age(
  cl_samples_mut_filtered, "mutational_burden", "age_prediction", "age_at_sampling", cancer_colors
)
