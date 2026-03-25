###################################
###   SUPPLEMENTARY_FIGURE_14   ###
#################################

# Generation of Supplementary Figure 14A and 14B

# ============================

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
library(readxl)

# set working directory
setwd("/group/iorio/Alessandro.D/epiclock")

# load scripts
source("src/global_params.R")
source("src/utils.R")
source("src/plot_generation_helpers.R")

# load enrichment data
KEGG_pan <- read_xlsx(paste0(results_path,
 "cell_lines/cell_lines_gene_dependencies/gsea_pancancer/KEGG_significant_pancancer_results.xlsx"))

# dotplot KEGG
GSEA_yo_dotplot(gsea_result = KEGG_pan, title = "Pan-Cancer Gene Enrichment KEGG",
                output_filename = "../../Supplementary Figures/Supplementary_Figure_14A.pdf")

# barplot KEGG
GSEA_yo_barplot(gsea_result = KEGG_pan, title = "Pan-Cancer Gene Enrichment KRGG",
                output_filename = "../../Supplementary Figures/Supplementary_Figure_14B.pdf")
