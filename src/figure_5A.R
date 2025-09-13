############################
###   FIGURE_5A         ###
###########################

# Generation of Figure 5A

# ============================
library(tidyverse)
library(dplyr)
library(data.table)
library(fgsea)
library(writexl)
library(MutExMatSorting)
library(reshape2)
library(readxl)

# set working directory
setwd("/group/iorio/Alessandro.D/epiclock")

# load scripts
source("src/global_params.R")
source("src/utils.R")
source("src/plot_generation_helpers.R")


# load significant target enrichments (BP) in different cancer types
significant_enrichments_bp_cancer_spec <- read_csv(paste0(results_path,
                                                "cell_lines/cell_lines_gene_dependencies/gsea_per_cancer_type/01_combined_BP_results.csv"))

# number of enrichments
for (cancer in names(table(significant_enrichments_bp_cancer_spec$cancer_type))) {
  print(paste0("There are ", table(significant_enrichments_bp_cancer_spec$cancer_type)[cancer],
               " gene target enrichments (GOBP) in ", cancer))
}

# dotplot for cancer specific GSEA (top 3 pathways per cancer type)
gsea_cancerspecific_dotplot(significant_enrichments_bp_cancer_spec,
                            title = NULL,
                            output_dir = "Figures/Figure5A.pdf")
