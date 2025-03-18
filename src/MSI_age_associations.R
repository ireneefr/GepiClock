####################################################
###   ANOVA MSI status                           ###
###   Author: Digilio Alessandro                 ###
###   Date: 05/03/2025                           ###
####################################################

# in this script is performed the correlation between the MSI status of the cell lines and
# the predicted age from the new model, and is compared with the results of the same analysis
# with donor age at sampling

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

# load data
# filter only samples with annotated MSI-MSS status
cl_samples_ms_status <- cl_samples_with_age_at_sampling %>%
  filter(msi_status %in% c("MSI", "MSS"))
print(paste("There are", nrow(cl_samples_ms_status),
            "samples with both 'age at sampling' and 'MSI status' correctly annotated"))

# consider only cell lines from cancer types with enough mss vs msi
cl_samples_MSIvsMSS <- cl_samples_ms_status[cl_samples_ms_status$cancer_type 
                                            %in% c("B-Lymphoblastic Leukemia",
                                                   "Colorectal Carcinoma",
                                                   "Gastric Carcinoma",
                                                   "Ovarian Carcinoma",
                                                   "T-Lymphoblastic Leukemia",
                                                   "Endometrial Carcinoma") ,]
table_MSI_status <- print_MSI_status(cl_samples_MSIvsMSS)

# ============================
# PAN-CANCER ANALYSIS
# ============================

# MSI STATUS AND AGE ASSOCIATION ANALYSIS (Pan-Cancer)
# -----------------------------------------------------
# we analyze whether cell lines classified as MSI (Microsatellite Instability) 
# or MSS (Microsatellite Stable) exhibit systematic differences in cell age 
# (both predicted age and donor age at sampling).

# compare MSI vs MSS using predicted age
test_result_pred_age_pan <- test_msi_age_association_pan_cancer(cl_samples_MSIvsMSS, "age_prediction")
boxplot_msi_age_pancancer(cl_samples_MSIvsMSS, "age_prediction",
                          p_value = test_result_pred_age_pan$p_value)

# compare MSI vs MSS using donor age at sampling
test_result_donor_age_pan <- test_msi_age_association_pan_cancer(cl_samples_MSIvsMSS, "age_at_sampling")
boxplot_msi_age_pancancer(cl_samples_MSIvsMSS, "age_at_sampling",
                          p_value = test_result_donor_age_pan$p_value)

# ============================
# CANCER-SPECIFIC ANALYSIS
# ============================

# MSI STATUS AND AGE ASSOCIATION ANALYSIS (Cancer Type Specific)
# -----------------------------------------------------
# we analyze whether cell lines classified as MSI (Microsatellite Instability) 
# or MSS (Microsatellite Stable) exhibit systematic differences in cell age 
# (both predicted age and donor age at sampling) among different cancer types (tested individually)

# compare MSI vs MSS using predicted age within cancer types
test_results_pred_age_cancer <- test_msi_age_association_cancer_specific(cl_samples_MSIvsMSS, "age_prediction")
boxplot_msi_age_cancer_specific(cl_samples_MSIvsMSS, "age_prediction")


# compare MSI vs MSS using donor age (age at sampling) within cancer types
test_results_donor_age_cancer <- test_msi_age_association_cancer_specific(cl_samples_MSIvsMSS, "age_at_sampling")
boxplot_msi_age_cancer_specific(cl_samples_MSIvsMSS, "age_at_sampling")

# landscape MSS/MSI plot (predicted age vs donor age)
# compute Cohen's d for predicted age
cohens_d_pred_age_cancer <- compute_cohens_d_cancer_specific(cl_samples_MSIvsMSS, "age_prediction")
# compute Cohen's d for donor age (age at sampling)
cohens_d_donor_age_cancer <- compute_cohens_d_cancer_specific(cl_samples_MSIvsMSS, "age_at_sampling")
# landscape plot
lollipop_cohens_d_cancer_specific(cohens_d_pred_age_cancer, cohens_d_donor_age_cancer)



