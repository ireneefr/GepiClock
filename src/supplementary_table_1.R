#################################
###   Supplementary Table 1   ###
#################################

source("src/utils.R")
source("src/global_params.R")
library(dplyr)

# Load TCGA samples
tcga_samples <- TCGA_samples(dir_data = "/Volumes/iorio/Irene/legacy/epiclock_old/data/")

# Get the TCGA project information
projects <- TCGAbiolinks::getGDCprojects()

# Merge TCGA samples with project info
tcga_samples_projects <- merge(tcga_samples, projects, by.x = "project", by.y = "id")

# Get Supplementary Table 1
df <- tcga_samples_projects %>%
  group_by(project, name) %>%  # Group by project
  summarise(number_samples = n(),
            min_age = min(age_at_index),
            max_age = max(age_at_index),
            median_age = median(age_at_index),
            .groups = "drop")  # Count rows
colnames(df) <- c("Project ID", "Project name", "Number of samples", "Min Age", "Max Age", "Median Age")
writexl::write_xlsx(df, "Supplementary Tables/Supplementary_Table_1.xlsx")
