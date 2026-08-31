##############################################
###   UMAP in TCGA samples                 ###
###   Author: Irene Fernández Rebollo      ###
###   Date: 24/08/2026                     ###
##############################################

#library(tidyverse)
library(TCGAbiolinks)
library(umap)
setwd("/group/iorio/Irene/epiclock_dev")
source("src/utils.R")

# Load data
tcga_bval <- TCGA_Bvalues()
tcga_samples <- TCGA_samples()
model_coefs <- read.csv("results/TCGA/model/model_coefs.csv")
cpgs_gepiclock <- model_coefs$X[model_coefs$s0 != 0][-1]

# CpGs variance
cpg_var <- apply(tcga_bval, 2, var, na.rm = TRUE)
top_cpgs <- names(sort(cpg_var, decreasing = TRUE))[1:50000]

# UMAP all CpGs
set.seed(42)
umap_result <- umap(tcga_bval[rownames(tcga_bval) %in% rownames(tcga_samples), top_cpgs])
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$project <- tcga_samples$project[match(rownames(umap_df), rownames(tcga_samples))]
umap_df$age_at_index <- tcga_samples$age_at_index[match(rownames(umap_df), rownames(tcga_samples))]

# UMAP GepiClock CpGs
set.seed(42)
umap_result_gepiclock <- umap(tcga_bval[rownames(tcga_bval) %in% rownames(tcga_samples), colnames(tcga_bval) %in% cpgs_gepiclock])
umap_df_gepiclock <- as.data.frame(umap_result_gepiclock$layout)
colnames(umap_df_gepiclock) <- c("UMAP1", "UMAP2")
umap_df_gepiclock$project <- tcga_samples$project[match(rownames(umap_df_gepiclock), rownames(tcga_samples))]
umap_df_gepiclock$age_at_index <- tcga_samples$age_at_index[match(rownames(umap_df_gepiclock), rownames(tcga_samples))]

# Save results
write.csv(umap_df, "results/umap/tcga_umap_coordinates.csv")
write.csv(umap_df_gepiclock, "results/umap/tcga_umap_gepiclock_coordinates.csv")
