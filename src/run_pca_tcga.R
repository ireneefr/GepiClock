#############################################
###   PCA in TCGA samples                 ###
###   Author: Irene Fernández Rebollo     ###
###   Date: 26/08/2026                    ###
#############################################

#library(tidyverse)
library(TCGAbiolinks)
setwd("/group/iorio/Irene/epiclock_dev")
source("src/utils.R")


# Load data
tcga_bval <- TCGA_Bvalues()
tcga_samples <- TCGA_samples()
model_coefs <- read.csv("results/TCGA/model/model_coefs.csv")
cpgs_gepiclock <- model_coefs$X[model_coefs$s0 != 0][-1]
tcga_palette <- read.csv("metadata/TCGAproject_palette.csv")
colnames(tcga_palette) <- c("project", "color")

# CpGs variance
cpg_var <- apply(tcga_bval, 2, var, na.rm = TRUE)
top_cpgs <- names(sort(cpg_var, decreasing = TRUE))[1:50000]


# PCA using the 50,000 most variable CpGs
pca_result <- prcomp(tcga_bval[,top_cpgs],
                     center = TRUE,
                     scale. = FALSE,
                     rank. = 2)

# PCA using GepiClock CpGs
pca_result_gepiclock <- prcomp(tcga_bval[,cpgs_gepiclock],
                               center = TRUE,
                               scale. = FALSE,
                               rank. = 2)

# Construct plotting data frames
pca_df <- as.data.frame(pca_result$x[, 1:2, drop = FALSE])
colnames(pca_df) <- c("PC1", "PC2")

pca_df$project <- tcga_samples[rownames(pca_df), "project"]
pca_df$age_at_index <- tcga_samples[rownames(pca_df), "age_at_index"]

pca_df_gepiclock <- as.data.frame(pca_result_gepiclock$x[, 1:2, drop = FALSE])
colnames(pca_df_gepiclock) <- c("PC1", "PC2")

pca_df_gepiclock$project <- tcga_samples[rownames(pca_df_gepiclock), "project"]

pca_df_gepiclock$age_at_index <- tcga_samples[rownames(pca_df_gepiclock), "age_at_index"]

# Percentage of variance explained
pca_variance <- 100 * pca_result$sdev^2 / sum(pca_result$sdev^2)
pca_variance_gepiclock <- 100 * pca_result_gepiclock$sdev^2 / sum(pca_result_gepiclock$sdev^2)

pc1_label <- sprintf("PC1 (%.1f%%)", pca_variance[1])
pc2_label <- sprintf("PC2 (%.1f%%)", pca_variance[2])

gepiclock_pc1_label <- sprintf("PC1 (%.1f%%)", pca_variance_gepiclock[1])
gepiclock_pc2_label <- sprintf("PC2 (%.1f%%)", pca_variance_gepiclock[2])

# Save coordinates
write.csv(pca_df, "results/pca/tcga_pca_coordinates.csv")
write.csv(pca_df_gepiclock, "results/pca/tcga_pca_gepiclock_coordinates.csv")

# Plot settings
cols <- setNames(tcga_palette$color, tcga_palette$project)

# PCA by TCGA project
pdf("results/pca/plots/pca_tcga_project.pdf", height = 8, width = 12)
ggplot(pca_df, aes(x = PC1, y = PC2, colour = project)) +
  geom_point() +
  scale_color_manual(values = cols, name = "Project") +
  labs(x = pc1_label, y = pc2_label) +
  theme_minimal(base_size = 20)
dev.off()

# GepiClock PCA by TCGA project
pdf("results/pca/plots/pca_tcga_gepiclock_project.pdf", height = 8, width = 12)
ggplot(pca_df_gepiclock, aes(x = PC1, y = PC2, colour = project)) +
    geom_point() +
    scale_color_manual(values = cols) +
    labs(x = gepiclock_pc1_label, y = gepiclock_pc2_label) +
    theme_minimal(base_size = 20)
dev.off()

# PCA by age
pdf("results/pca/plots/pca_tcga_age.pdf", height = 8, width = 10)
ggplot(pca_df, aes(x = PC1, y = PC2, colour = age_at_index)) +
    geom_point() +
    scale_color_viridis_c(name = "Age", na.value = "grey80") +
    labs(x = pc1_label, y = pc2_label) +
    theme_minimal(base_size = 20)
dev.off()

# GepiClock PCA by age
pdf("results/pca/plots/pca_tcga_gepiclock_age.pdf", height = 8, width = 10)
ggplot(pca_df_gepiclock, aes(x = PC1, y = PC2, colour = age_at_index)) +
    geom_point() +
    scale_color_viridis_c(name = "Age", na.value = "grey80") +
    labs(x = gepiclock_pc1_label, y = gepiclock_pc2_label) +
    theme_minimal(base_size = 20)
dev.off()