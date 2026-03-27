##################################
###   Supplementary Figure 3   ###
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
tcga_samples_projects$sample_type <- factor(tcga_samples_projects$sample_type,
                                            levels = c("Solid Tissue Normal",
                                                       "Primary Tumor",
                                                       "Primary Blood Derived Cancer - Peripheral Blood"))

# Get Supplementary Figure 3
ggplot(tcga_samples_projects, aes(x = reorder(project, project, length, decreasing = TRUE), fill = sample_type)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 3) +
  scale_fill_manual(values = c("chartreuse3", "steelblue1", "steelblue3")) +
  xlab("") + ylab("Number of samples") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = "white", color = "white", size = 0.3)) +
  guides(fill = guide_legend(title = "Sample type"))
ggsave("Supplementary Figures/Supplementary_Figure_3.pdf")
