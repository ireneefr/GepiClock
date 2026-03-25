###################################
###   SUPPLEMENTARY_FIGURE_12   ###
#################################

# Generation of Supplementary Figure 12A and 12B

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
GO_CC_pan <- read_xlsx(paste0(results_path,
 "cell_lines/cell_lines_gene_dependencies/gsea_pancancer/GOCC_significant_pancancer_results.xlsx"))

# dotplot GO Cellular Component (CC)
GSEA_yo_dotplot(gsea_result = GO_CC_pan, title = "Pan-Cancer Gene Enrichment CC",
                output_filename = "../../Supplementary Figures/Supplementary_Figure_12A.pdf")

# barplot GO Cellular Component (CC)
GSEA_yo_barplot(gsea_result = GO_CC_pan, title = "Pan-Cancer Gene Enrichment CC",
                output_filename = "../../Supplementary Figures/Supplementary_Figure_12B.pdf")
