####################################################
###   Data Availability                          ###
###   Author: Digilio Alessandro                 ###
###   Date: 05/03/2025                           ###
####################################################

# in this script are collected all the objects containing data and metadata

source("src/utils.R")

library(tidyverse)
library(readxl)

# ==================================================
# MODEL'S CpGs & COEFFICIENTS
# ==================================================

# load coefficients/CpGs of the new model
model <- load_model_coefficients("results/model/model_coefs.csv")
model_CpGs <- model$model_CpGs
coefficients_new_model <- model$coefficients_new_model
dim(coefficients_new_model)  
length(model_CpGs)

# ==================================================
# CELL LINES' METHYLATION DATA & ANNOTATIONS
# ==================================================

# load the model methylation data (b-values) and the annotations from CMP of the cell lines
cell_lines_data <- load_b_values_and_annotations(
  b_values_path = "data/b_values/CLs_methylation_data.csv", # load methylation data
  model_CpGs = model_CpGs,# load model's CpGs
  iorio_path = "metadata/CMP_annotations/IorioCell2016-MethylAccesionCellLines.txt", # load iorio cls ACcession code 
  annotations_path = "metadata/CMP_annotations/model_list_20240110.csv", # load cls metadata (annotations)
  experimental_settings_path = "metadata/CMP_GROWTH_20250114.txt" # load doubling time annotations (experimental settings)
   )
# extract matched data
b_values <- cell_lines_data$b_values
cl_samples <- cell_lines_data$cl_samples
# dimensions
dim(b_values)  # 4862 CpGs, 920 cell lines
print(head(b_values[, 1:6]))
dim(cl_samples)  # 920 cell lines
print(head(cl_samples[, 1:6]))
# write.csv(b_values, "data/b_values/new_model_b_values.csv", row.names = TRUE))

# ==================================================
# PREDICTED AGE ANNOTATIONS
# ==================================================

# add predicted age in the metadata
source("src/model_age_prediction_CLs.R")

# ==================================================
# AGE AT SAMPLING ANNOTATIONS
# ==================================================

# filter only cl_samples with age at sampling annotations
cl_samples_with_age_at_sampling <- cl_samples[!is.na(cl_samples$age_at_sampling) , ] 
print((paste("The total number of samples is", nrow(cl_samples), ",while the number of samples with annotated 'age at sampling' is",
             nrow(cl_samples_with_age_at_sampling))))

# ==================================================
# CELL LINES' DRUG RESPONSE DATA (GDSC)
# ==================================================

# drug response data from the GDSC databank
GDSC1_IC50s <- read_excel("data/GDSC/GDSC1_fitted_dose_response_27Oct23.xlsx")
GDSC2_IC50s <- read_excel("data/GDSC/GDSC2_fitted_dose_response_27Oct23.xlsx")
# combine and clean GDSC Data
GDSC_combined <- bind_rows(GDSC1_IC50s, GDSC2_IC50s) %>%
  distinct(COSMIC_ID, DRUG_NAME, CELL_LINE_NAME, .keep_all = TRUE) # distinct: duplicate rows are removed
print("The GDSC1 and GDSC2 datasets have been loaded and merged")

# ==================================================
# CELL LINES' GENE DEPENDENCY DATA (DepMap)
# ==================================================

# gene effect data from DepMap (pan-cancer)
gene_effect_and_cl_filtered <- load_and_filter_genedep_pancancer(
  depmap_file = "data/depmap/CRISPRGeneEffect.csv", # Gene effect matrix from DepMap Public 24Q2 Files
  essential_genes_file = "data/depmap/AchillesCommonEssentialControls.csv", # Achilles Common Essential genes
  nonessential_genes_file = "data/depmap/AchillesNonessentialControls.csv", # Achilles Never Essential genes
  bagel_data_path = "data/BAGEL/", # additional essential & non-essential gene sets (BAGEL dataset)
  cl_samples = cl_samples, # cell line metadata for filtering
  min_cell_lines = 2 # retain genes depleted (-0.5) in at least 2 cell lines
)
# extract data
gene_effect_filtered <- gene_effect_and_cl_filtered$gene_effect_filtered
cl_samples_gene_filtered <- gene_effect_and_cl_filtered$cl_samples_gene_filtered
print(dim(gene_effect_filtered)) 
print(dim(cl_samples_gene_filtered))

# ==================================================
# CANCER TYPE-SPECIFIC GENE DEPENDENCY DATA (DepMap)
# ==================================================

# gene effect data from DepMap (cancer-type specific)
gene_effect_and_cl_filtered_cancer_specific <- load_and_filter_genedep_cancer_type(
  depmap_file = "data/depmap/CRISPRGeneEffect.csv", # Gene effect matrix from DepMap Public 24Q2 Files
  essential_genes_file = "data/depmap/AchillesCommonEssentialControls.csv", # Achilles Common Essential genes
  nonessential_genes_file = "data/depmap/AchillesNonessentialControls.csv", # Achilles Never Essential genes
  bagel_data_path = "data/BAGEL/", # additional essential & non-essential gene sets (BAGEL dataset)
  cl_samples = cl_samples,  # cell line metadata for filtering
  min_cell_lines = 2  # retain genes depleted (-0.5) in at least 2 cell lines
)
# check the number of genes for each cancer type
gene_counts_per_cancer <- sapply(gene_effect_and_cl_filtered_cancer_specific, nrow)
print(gene_counts_per_cancer)
print("All the DepMap dependency data have been loaded, filtered and scaled")






