#################################
###   Supplementary Table 1   ###
#################################

source("src/utils.R")
source("src/global_params.R")
library(dplyr)

# Load TCGA samples
tcga_samples <- TCGA_samples(dir_data = "/Volumes/iorio/Irene/legacy/epiclock_old/data/")
dim(tcga_samples)
tcga_samples[1:5,1:5]
colnames(tcga_samples)

# Get the TCGA project information
projects <- TCGAbiolinks::getGDCprojects()
head(projects)

# Merge TCGA samples with project info
tcga_samples_projects <- merge(tcga_samples, projects, by.x = "project", by.y = "id")

# Get Supplementary Table 1
df <- tcga_samples_projects %>%
  group_by(project, name) %>%  # Group by project
  summarise(number_samples = n(), .groups = "drop")  # Count rows
colnames(df) <- c("Project ID", "Project name", "Number of samples")
writexl::write_xlsx(df, "Supplementary Tables/Supplementary_Table_1.xlsx")
