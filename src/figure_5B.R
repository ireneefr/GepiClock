############################
###   FIGURE_5B         ###
###########################

# Generation of Figure 5B_1 and 5C_specific_genes

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
source("src/load_data.R")

# load DepMap x cell lines cancer specific data (filtered after category selection)
gene_dep_age_filtered <- readRDS(paste0(results_path, "cell_lines/cell_lines_gene_dependencies/cell_lines_per_genes_cancer_specific.rds"))

# perform anova anova model for gene dependencies
cancer_specific_gene_anova_results <- anova_genedep_cancer_specific(
  scaled_DepMat_list = gene_effect_and_cl_filtered_cancer_specific, 
  cl_samples_filtered = cl_samples_gene_filtered
)
summarize_significant_genes(cancer_specific_gene_anova_results)

# extract and plot cumulative correlation results across cancer types
cumulative_correlation_results_geneDep <- read_xlsx(paste0(results_path,
                                                           "cell_lines/cell_lines_gene_dependencies/anova_cancer_specific_cumulative_correlation_results.xlsx"))

# identify cancer types with statistically significant age-gene correlation
# - Acute Myeloid Leukemia
# - Bladder Carcinoma 
# - Other Solid Cancers 

# extract list of significant genes per cancer type
significant_genes_list <- extract_significant_genes(cancer_specific_gene_anova_results)
# plot volcano plot for selected cancer types
plot_cancer_specific_geneDep_volcano(cumulative_correlation_results_geneDep, 
                                     c("Acute Myeloid Leukemia", 
                                       "Bladder Carcinoma", 
                                       "Other Solid Cancers"),
                                     output_dir = "Figures/Figure5B1.pdf" )

# individual correlation plots for statistically significant genes
plot_significant_genes(gene_dep_age_filtered, significant_genes_list,
                       output_dir = "Figures/")

