###################################
###   SUPPLEMENTARY_FIGURE_10   ###
#################################

# Generation of Supplementary Figure 10

# ============================
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
print(paste("In the pan-cancer analysis, the number of cell lines after the category selection is",
            length(unique(gene_dep_age$BROAD_ID))))

# plot sample categories
depmap_combined <- merged_dep_cls$depmap_combined
plot_sample_categories_pan_gene(depmap_combined, thresh_number = 10,
                                output_dir = "Supplementary Figures/Supplementary_Figure_10.pdf")
