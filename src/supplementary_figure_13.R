###################################
###   SUPPLEMENTARY_FIGURE_13   ###
#################################

# Generation of Supplementary Figure 13A and 13B

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
GO_MF_pan <- read_xlsx(paste0(results_path,
 "cell_lines/cell_lines_gene_dependencies/gsea_pancancer/GOMF_significant_pancancer_results.xlsx"))

# dotplot GO Molecular Function (MF)
GSEA_yo_dotplot(gsea_result = GO_MF_pan, title = "Pan-Cancer Gene Enrichment MF",
                output_filename = "../../Supplementary Figures/Supplementary_Figure_13A.pdf")

# barplot GO Molecular Function (MF)
GSEA_yo_barplot(gsea_result = GO_MF_pan, title = "Pan-Cancer Gene Enrichment MF",
                output_filename = "../../Supplementary Figures/Supplementary_Figure_13B.pdf")
