####################################################
###   ANOVA genetic dependencies                 ###
###   Author: Digilio Alessandro                 ###
###   Date: 05/03/2025                           ###
####################################################

# in this script is performed the correlation between the gene dependencies of the cell lines 
# in terms of "gene effect" (essentality given from a CRISPR screen) data from depmap (corrected with ANOVA
# from the variability of tissue/cancer type and  MSI status) and the predicted age from the new model

library(tidyverse)
library(dplyr)
library(data.table)
library(fgsea)
library(writexl)
library(GO.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)

# set working directory
setwd("/group/iorio/Alessandro.D/epiclock")

# load scripts
source("src/global_params.R")
source("src/utils.R")
source("src/plot_generation_helpers.R")
source("src/load_data.R")

# set seed for reproducibility
set.seed(123)

# ============================
# DATA LOADING & PREPROCESSING
# ============================
# generate merged dataset for gene dependencies
merged_dep_cls <- generate_merged_gene_dataset(gene_effect_filtered, cl_samples_gene_filtered, threshold = 10)
gene_dep_age <- merged_dep_cls$gene_dep_age
print(head(gene_dep_age))
print(paste("In the pan-cancer analysis, the number of cell lines after the category selection is", length(unique(gene_dep_age$BROAD_ID))))

# save data
saveRDS(gene_dep_age, "results/cell_lines/cell_lines_gene_dependencies/cell_lines_per_genes_pancancer.rds")

# plot sample categories
depmap_combined <- merged_dep_cls$depmap_combined
plot_sample_categories_pan_gene(depmap_combined, thresh_number = 10,
                                output_dir = paste0(figures_path, "tissue_or_cancer_selection_gene_pan.png"))
plot_included_tissue_pan_gene(depmap_combined)

# ============================
# PAN-CANCER ANALYSIS
# ============================

# ANOVA MODEL FOR GENETIC DEPENDENCIES (Pan-Cancer)
# ------------------------------------------
# here is performed an ANOVA analysis to evaluate how tissue/cancer type and MSI status influence gene dependencies (gene effect) across cell lines.
# The model accounts for the effects of these factors and extracts residuals.Residuals represent the unexplained variance in gene dependency after 
# considering tissue/cancer type and MSI status.A correlation analysis is performed between residuals and predicted age to assess whether age 
# influences gene dependencies.False Discovery Rate (FDR) correction is applied to control for multiple comparisons, identifying genes with significant 
# age-response associations.

# perform anova anova model for gene dependencies
anova_gene_pancancer_results <- anova_analysis_genedep_pan(gene_dep_age)
# plot correlation between predicted age and residual gene sensitivity
correlation_results_ANOVA_gene_pan <- anova_gene_pancancer_results$correlation_results

# save results
write_xlsx(correlation_results_ANOVA_gene_pan, "results/cell_lines/cell_lines_drug_sensitivity/anova_pancancer_correlation_results.xlsx")

plot_predage_gene_effect_correlation_pan_gene(correlation_results_ANOVA_gene_pan)

# GENE SET ENRICHMENT ANALYSIS (GSEA) - PAN-CANCER
# ------------------------------------------------
# perform Gene Set Enrichment Analysis (GSEA) using genes ranked by their correlation 
# between predicted age and residual gene dependency (gene_effect).
# identifies enriched pathways associated with age-dependent gene dependencies
gsea_results_pan_cancer <- run_gsea_pan_cancer(correlation_results_ANOVA_gene_pan)

# significant enrichment
GO_BP_pan <- gsea_results_pan_cancer$GO_BP
print(GO_BP_pan)
GO_MF_pan <- gsea_results_pan_cancer$GO_MF
print(GO_MF_pan)
GO_CC_pan <- gsea_results_pan_cancer$GO_CC
print(GO_CC_pan)
KEGG_pan <- gsea_results_pan_cancer$KEGG
print(KEGG_pan)

# dot plots
# GO Biological Process (BP)
GSEA_yo_dotplot(gsea_result = GO_BP_pan, title = "Pathway Enrichment by Predicted Group", output_filename = "gsea_gobp_pan.png")
# GO Molecular Function (MF)
GSEA_yo_dotplot(gsea_result = GO_MF_pan, title = "Pathway Enrichment by Predicted Group", output_filename = "gsea_gomf_pan.png")
# GO Cellular Component (CC)
GSEA_yo_dotplot(gsea_result = GO_CC_pan, title = "Pathway Enrichment by Predicted Group", output_filename = "gsea_gocc_pan.png")
# KEGG Pathways
GSEA_yo_dotplot(gsea_result = KEGG_pan, title = "Pathway Enrichment by Predicted Group", output_filename = "gsea_kegg_pan.png")

# ============================
# CANCER TYPE-SPECIFIC ANALYSIS
# ============================

# filter dataset to include only cell lines categorized as cancer types with sufficient samples (≥10)
gene_dep_age_filtered <- gene_dep_age[gene_dep_age$cancer_specific == "Included as Cancer Type", ]
gene_dep_age_filtered$residuals <- NULL
head(gene_dep_age_filtered)
print(paste("In the cancer type-specific analysis, the number of cell lines after the category selection is", length(unique(gene_dep_age_filtered$BROAD_ID))))

# save data
saveRDS(gene_dep_age_filtered, "results/cell_lines/cell_lines_gene_dependencies/cell_lines_per_genes_cancer_specific.rds")

# plot sample categories
plot_sample_categories_cancer_specific_gene(gene_dep_age_filtered,
                                            output_dir = paste0(figures_path,
                                                                "cancer_selection_genedep_canc_spec.png")

# ANOVA MODEL FOR GENE DEPENDENCIES (Cancer-Specific)
# -----------------------------------------------
# performs the same ANOVA model within each cancer type to analyze gene dependencies, extracts residuals
# and assesses correlation with predicted age.
# the aim is the identification of cancer types where gene dependencies is significantly age-associated.

# perform anova anova model for gene dependencies
cancer_specific_gene_anova_results <- anova_genedep_cancer_specific(
  scaled_DepMat_list = gene_effect_and_cl_filtered_cancer_specific, 
  cl_samples_filtered = cl_samples_gene_filtered
)
summarize_significant_genes(cancer_specific_gene_anova_results)

# extract and plot cumulative correlation results across cancer types
cumulative_correlation_results_geneDep <- cancer_specific_gene_anova_results$cumulative_correlation_results

# save results
write_xlsx(cumulative_correlation_results_geneDep, "results/cell_lines/cell_lines_gene_dependencies/anova_cancer_specific_cumulative_correlation_results.xlsx")

plot_predage_geneDep_cumulative_correlation_cancer_specific(cumulative_correlation_results_geneDep)

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
                                     output_dir = paste0(figures_path, "volcano_cancer_specific_geneDep.png"))

# individual correlation plots for statistically significant genes
plot_significant_genes(gene_dep_age_filtered, significant_genes_list, output_dir = paste0(results_path, "genes/"))

# GENE SET ENRICHMENT ANALYSIS (GSEA) - CANCER-SPECIFIC
# -----------------------------------------------------
# performs GSEA for each cancer type separately.
# ranks genes based on correlation between predicted age and gene dependency (gene_effect).
# identifies enriched gene targets specific to individual cancer types.
# run GSEA for gene dependencies in cancer-specific analysis
# this step can take a while...
gsea_results_cancer_specific <- run_gsea_cancer_specific(cumulative_correlation_results_geneDep)

# extract significant target enrichments (BP) in different cancer types
significant_enrichments_bp_cancer_spec <- gsea_results_cancer_specific$BP
# number of enrichments
for (cancer in names(table(significant_enrichments_bp_cancer_spec$cancer_type))) {
  print(paste0("There are ", table(significant_enrichments_bp_cancer_spec$cancer_type)[cancer],
               " gene target enrichments (GOBP) in ", cancer))
}

# dotplot for cancer specific GSEA (top 3 pathways per cancer type)
gsea_cancerspecific_dotplot(significant_results_bp_df, title = NULL, output_dir = paste0(figures_path, "gsea_dotplot_cancer_specific_MutExMatSorting.png"))









