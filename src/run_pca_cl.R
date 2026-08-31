#############################################
###   PCA in Cell Lines samples           ###
###   Author: Irene Fernández Rebollo     ###
###   Date: 26/08/2026                    ###
#############################################

#library(tidyverse)
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


# PCA using the 50,000 most variable CpGs
pca_result <- prcomp(cl_bval[,top_cpgs],
                     center = TRUE,
                     scale. = FALSE,
                     rank. = 2)

# PCA using GepiClock CpGs
pca_result_gepiclock <- prcomp(cl_bval[,cpgs_gepiclock],
                               center = TRUE,
                               scale. = FALSE,
                               rank. = 2)

# Construct plotting data frames
pca_df <- as.data.frame(pca_result$x[, 1:2, drop = FALSE])
colnames(pca_df) <- c("PC1", "PC2")

pca_df$project <- cl_samples$tissue[match(rownames(pca_df), cl_samples$CAccession)]
pca_df$age_at_index <- cl_samples$age_at_sampling[match(rownames(pca_df), cl_samples$CAccession)]

pca_df_gepiclock <- as.data.frame(pca_result_gepiclock$x[, 1:2, drop = FALSE])
colnames(pca_df_gepiclock) <- c("PC1", "PC2")

pca_df_gepiclock$project <- cl_samples$tissue[match(rownames(pca_df_gepiclock), cl_samples$CAccession)]
pca_df_gepiclock$age_at_index <- cl_samples$age_at_sampling[match(rownames(pca_df_gepiclock), cl_samples$CAccession)]

# Percentage of variance explained
pca_variance <- 100 * pca_result$sdev^2 / sum(pca_result$sdev^2)
pca_variance_gepiclock <- 100 * pca_result_gepiclock$sdev^2 / sum(pca_result_gepiclock$sdev^2)

pc1_label <- sprintf("PC1 (%.1f%%)", pca_variance[1])
pc2_label <- sprintf("PC2 (%.1f%%)", pca_variance[2])

gepiclock_pc1_label <- sprintf("PC1 (%.1f%%)", pca_variance_gepiclock[1])
gepiclock_pc2_label <- sprintf("PC2 (%.1f%%)", pca_variance_gepiclock[2])

# Save coordinates
write.csv(pca_df, "results/pca/cl_pca_coordinates.csv")
write.csv(pca_df_gepiclock, "results/pca/cl_pca_gepiclock_coordinates.csv")

# Plot settings
cols <- tissue_colors

# PCA by CL tissue
pdf("results/pca/plots/pca_cl_tissue.pdf", height = 8, width = 12)
ggplot(pca_df, aes(x = PC1, y = PC2, colour = project)) +
  geom_point(size = 3) +
  scale_color_manual(values = cols, name = "Tissue") +
  labs(x = pc1_label, y = pc2_label) +
  theme_minimal(base_size = 20)
dev.off()

# GepiClock PCA by CL tissue
pdf("results/pca/plots/pca_cl_gepiclock_tissue.pdf", height = 8, width = 12)
ggplot(pca_df_gepiclock, aes(x = PC1, y = PC2, colour = project)) +
  geom_point(size = 3) +
  scale_color_manual(values = cols, name = "Tissue") +
  labs(x = gepiclock_pc1_label, y = gepiclock_pc2_label) +
  theme_minimal(base_size = 20)
dev.off()

# PCA by age
pdf("results/pca/plots/pca_cl_age.pdf", height = 8, width = 10)
ggplot(pca_df, aes(x = PC1, y = PC2, colour = age_at_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(name = "Age", na.value = "grey80") +
  labs(x = pc1_label, y = pc2_label) +
  theme_minimal(base_size = 20)
dev.off()

# GepiClock PCA by age
pdf("results/pca/plots/pca_cl_gepiclock_age.pdf", height = 8, width = 10)
ggplot(pca_df_gepiclock, aes(x = PC1, y = PC2, colour = age_at_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(name = "Age", na.value = "grey80") +
  labs(x = gepiclock_pc1_label, y = gepiclock_pc2_label) +
  theme_minimal(base_size = 20)
dev.off()