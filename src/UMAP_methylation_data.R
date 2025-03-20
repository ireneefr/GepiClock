####################################################
###   UMAP methylation data                      ###
###   Author: Digilio Alessandro                 ###
###   Date: 05/03/2025                           ###
####################################################

# this script performs UMAP analysis on methylation data (b-values) of cell lines 
# to visualize clustering patterns based on tissue of origin, cancer type, and MSI status.

library(ggplot2)
library(umap)
library(tidyverse)
library(data.table)
library(readxl)

# set working directory
setwd("/group/iorio/Alessandro.D/epiclock")

# load scripts
source("src/global_params.R")
source("src/utils.R")
source("src/plot_generation_helpers.R")

# ==================================================
# LOAD & PROCESS METHYLATION DATA (450k CpGs)
# ==================================================

# load 450k methylation data (entire methylation profile of the cell lines)
all_methylation_data <- load_450k_methylation_data(
  all_b_values_path = paste0(data_path, "b_values/CLs_methylation_data.csv"),  # complete methylation data
  iorio_path = paste0(metadata_path, "CMP_annotations/IorioCell2016-MethylAccesionCellLines.txt"), # Iorio dataset
  annotations_path = paste0(metadata_path, "CMP_annotations/model_list_20240110.csv"), # cell line metadata
  experimental_settings_path = paste0(metadata_path, "CMP_GROWTH_20250114.txt") # doubling time annotations
)

# extract processed datasets
all_b_values <- all_methylation_data$all_b_values   # CpG methylation data 
cl_samples <- all_methylation_data$cl_samples       # cell line metadata 

# check data dimensions
dim(all_b_values)   
dim(cl_samples)

# ==================================================
# UMAP ANALYSIS ON COMPLETE METHYLOME
# ==================================================

# define output paths
output_tissue <- paste0(figures_path, "umap/umap_all_b_values_tissue.png")
output_cancer <- paste0(figures_path, "umap/umap_all_b_values_cancer_type.png")
output_msi <- paste0(figures_path, "umap/umap_all_b_values_msi.png")

# generate UMAP plots
plot_umap_methylation_data(all_b_values, cl_samples, "tissue", output_tissue)       # by tissue of origin
plot_umap_methylation_data(all_b_values, cl_samples, "cancer_type", output_cancer)  # by cancer type
plot_umap_methylation_data(all_b_values, cl_samples, "msi_status", output_msi)      # by MSI status

# ==================================================
# LOAD & PROCESS MODEL METHYLATION DATA (4862 CpGs)
# ==================================================

# define output paths
output_tissue_model <- paste0(figures_path, "umap/umap_b_values_tissue.png")
output_cancer_model <- paste0(figures_path, "umap/umap_b_values_cancer_type.png")
output_msi_model <- paste0(figures_path, "umap/umap_b_values_msi.png")

# load coefficients/CpGs of the new model
model <- load_model_coefficients("results/model/model_coefs.csv")
model_CpGs <- model$model_CpGs
coefficients_new_model <- model$coefficients_new_model
dim(coefficients_new_model)  
length(model_CpGs)

# load model-specific methylation data (b-values) & annotations
model_methylation_data <- load_b_values_and_annotations(
  b_values_path = paste0(data_path, "b_values/CLs_methylation_data.csv"),  # methylation data
  model_CpGs = model_CpGs,  # selected CpGs from the model
  iorio_path = paste0(metadata_path, "CMP_annotations/IorioCell2016-MethylAccesionCellLines.txt"), # Iorio dataset
  annotations_path = paste0(metadata_path, "CMP_annotations/model_list_20240110.csv"), # cell line metadata
  experimental_settings_path = paste0(metadata_path, "CMP_GROWTH_20250114.txt") # doubling time annotations
)

# extract processed datasets
b_values <- model_methylation_data$b_values  # filtered methylation data
cl_samples <- model_methylation_data$cl_samples  # matched metadata

# ==================================================
# UMAP ANALYSIS ON MODEL-SPECIFIC CpGs
# ==================================================

# generate UMAP plots
plot_umap_methylation_data(b_values, cl_samples, "tissue", output_tissue_model)       # by tissue of origin
plot_umap_methylation_data(b_values, cl_samples, "cancer_type", output_cancer_model)  # by cancer type
plot_umap_methylation_data(b_values, cl_samples, "msi_status", output_msi_model)      # by MSI status
