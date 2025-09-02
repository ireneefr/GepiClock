##################################
###   Supplementary Figure 1   ###
##################################

source("src/utils.R")
source("src/global_params.R")
library(ggplot2)

# Load TCGA samples
tcga_samples <- TCGA_samples(dir_data = "/Volumes/iorio/Irene/legacy/epiclock_old/data/")

# Get the TCGA project information
projects <- TCGAbiolinks::getGDCprojects()

# Merge TCGA samples with project info
tcga_samples_projects <- merge(tcga_samples, projects, by.x = "project", by.y = "id")

# Get Supplementary Figure 1
tcga_palette <- read.csv("metadata/TCGAproject_palette.csv")
ggplot(tcga_samples_projects, aes(x = age_at_index, y = project, fill = project)) +
  geom_boxplot() +
  scale_fill_manual(values = setNames(tcga_palette[,2], tcga_palette[,1])) +
  xlab("Age") + ylab("") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none")
ggsave("Supplementary Figures/Supplementary_Figure_1.png")
