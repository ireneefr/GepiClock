##############################################
###   UMAP in Cell Line samples            ###
###   Author: Irene Fernández Rebollo      ###
###   Date: 24/08/2026                     ###
##############################################

#library(tidyverse)
library(TCGAbiolinks)
library(umap)
setwd("/group/iorio/Irene/epiclock_dev")
source("src/utils.R")
source("src/global_params.R")


# load 450k methylation data (entire methylation profile of the cell lines)
all_methylation_data <- load_450k_methylation_data(
  all_b_values_path = paste0(data_path, "b_values/CLs_methylation_data.csv"),  # complete methylation data
  iorio_path = paste0(metadata_path, "CMP_annotations/IorioCell2016-MethylAccesionCellLines.txt"), # Iorio dataset
  annotations_path = paste0(metadata_path, "CMP_annotations/model_list_20240110.csv"), # cell line metadata
  experimental_settings_path = paste0(metadata_path, "CMP_GROWTH_20250114.txt") # doubling time annotations
)
cl_bval <- t(all_methylation_data$all_b_values)
cl_samples <- all_methylation_data$cl_samples
model_coefs <- read.csv("results/TCGA/model/model_coefs.csv")
cpgs_gepiclock <- model_coefs$X[model_coefs$s0 != 0][-1]

# CpGs variance
cpg_var <- apply(cl_bval, 2, var, na.rm = TRUE)
top_cpgs <- names(sort(cpg_var, decreasing = TRUE))[1:50000]

# UMAP all CpGs
set.seed(42)
umap_result <- umap(cl_bval[rownames(cl_bval) %in% cl_samples$CAccession, top_cpgs])
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$project <- cl_samples$tissue[match(rownames(umap_df), cl_samples$CAccession)]
umap_df$age_at_index <- cl_samples$age_at_sampling[match(rownames(umap_df), cl_samples$CAccession)]

# UMAP GepiClock CpGs
set.seed(42)
umap_result_gepiclock <- umap(cl_bval[rownames(cl_bval) %in% cl_samples$CAccession, colnames(cl_bval) %in% cpgs_gepiclock])
umap_df_gepiclock <- as.data.frame(umap_result_gepiclock$layout)
colnames(umap_df_gepiclock) <- c("UMAP1", "UMAP2")
umap_df_gepiclock$project <- cl_samples$tissue[match(rownames(umap_df_gepiclock), cl_samples$CAccession)]
umap_df_gepiclock$age_at_index <- cl_samples$age_at_sampling[match(rownames(umap_df_gepiclock), cl_samples$CAccession)]

# Save results
write.csv(umap_df, "results/umap/cl_umap_coordinates.csv")
write.csv(umap_df_gepiclock, "results/umap/cl_umap_gepiclock_coordinates.csv")
