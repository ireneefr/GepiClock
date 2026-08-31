##############################################
###   Plot UMAP in TCGA and Cell Lines     ###
###   Author: Irene Fernández Rebollo      ###
###   Date: 24/08/2026                     ###
##############################################

#library(tidyverse)
library(ggplot2)
setwd("/group/iorio/Irene/epiclock_dev")
source("src/global_params.R")
tcga_palette <- read.csv("metadata/TCGAproject_palette.csv")
colnames(tcga_palette) <- c("project", "color")
umap_df <- read.csv("results/umap/tcga_umap_coordinates.csv")
umap_df_gepiclock <- read.csv("results/umap/tcga_umap_gepiclock_coordinates.csv")
cl_umap_df <- read.csv("results/umap/cl_umap_coordinates.csv")
cl_umap_df_gepiclock <- read.csv("results/umap/cl_umap_gepiclock_coordinates.csv")

# Plot UMAP TCGA
cols <- setNames(tcga_palette$color, tcga_palette$project)
pdf("results/umap/plots/umap_tcga_project.pdf", height = 8, width = 12)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = project)) +
  geom_point() +
  scale_color_manual(values = cols) +
  theme_minimal(base_size = 20)
dev.off()
pdf("results/umap/plots/umap_tcga_gepiclock_project.pdf", height = 8, width = 12)
ggplot(umap_df_gepiclock, aes(x = UMAP1, y = UMAP2, colour = project)) +
  geom_point() +
  scale_color_manual(values = cols) +
  theme_minimal(base_size = 20)
dev.off()

pdf("results/umap/plots/umap_tcga_age.pdf", height = 8, width = 10)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = age_at_index)) +
  geom_point() +
  scale_color_viridis_c(name = "Age") +
  theme_minimal(base_size = 20)
dev.off()
pdf("results/umap/plots/umap_tcga_gepiclock_age.pdf", height = 8, width = 10)
ggplot(umap_df_gepiclock, aes(x = UMAP1, y = UMAP2, colour = age_at_index)) +
  geom_point() +
  scale_color_viridis_c(name = "Age") +
  theme_minimal(base_size = 20)
dev.off()

# Plot UMAP CL
cols <- tissue_colors
pdf("results/umap/plots/umap_cl_tissue.pdf", height = 8, width = 14)
ggplot(cl_umap_df, aes(x = UMAP1, y = UMAP2, colour = project)) +
  geom_point() +
  scale_color_manual(values = cols, name = "Tissue") +
  theme_minimal(base_size = 20)
dev.off()
pdf("results/umap/plots/umap_cl_gepiclock_tissue.pdf", height = 8, width = 14)
ggplot(cl_umap_df_gepiclock, aes(x = UMAP1, y = UMAP2, colour = project)) +
  geom_point() +
  scale_color_manual(values = cols, name = "Tissue") +
  theme_minimal(base_size = 20)
dev.off()

pdf("results/umap/plots/umap_cl_age.pdf", height = 8, width = 10)
ggplot(cl_umap_df, aes(x = UMAP1, y = UMAP2, colour = age_at_index)) +
  geom_point() +
  scale_color_viridis_c(name = "Age") +
  theme_minimal(base_size = 20)
dev.off()
pdf("results/umap/plots/umap_cl_gepiclock_age.pdf", height = 8, width = 10)
ggplot(cl_umap_df_gepiclock, aes(x = UMAP1, y = UMAP2, colour = age_at_index)) +
  geom_point() +
  scale_color_viridis_c(name = "Age") +
  theme_minimal(base_size = 20)
dev.off()
