############################
###   FIGURE_4B         ###
###########################

# Generation of Figure 4B and 4B2

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

# load GDSC_age object (drug sensitivity + putative targets)
GDSC_age_res <- readRDS(paste0(results_path, "cell_lines/cell_lines_drug_sensitivity/cell_lines_per_drugs_cancer_specific.rds"))

# check
print(colnames(GDSC_age_res)[1:20])
stopifnot(all(c("DRUG_NAME", "PUTATIVE_TARGET") %in% colnames(GDSC_age_res)))

# load enrichment data
dsea_cancer_specific_results <- read_xlsx(paste0(results_path,
                                            "cell_lines/cell_lines_drug_sensitivity/dsea_per_cancer_type/significant_enrichment_dsea_cancerspecific.xlsx"))

# load correlation results used to rank drugs
correlation_results_ANOVA <- read_xlsx(
  paste0(results_path, "cell_lines/cell_lines_drug_sensitivity/anova_cancer_specific_cumulative_correlation_results.xlsx")
)
print(head(correlation_results_ANOVA))

# plot figure 4B (drug target enrichment - cancer specific)
# plot cancer-type-specific DSEA results
dsea_cancerspecific_dotplot(dsea_cancer_specific_results, "Cancer-Specific DSEA",
                            output_dir = "Figures/Figure4B.pdf")

# plot figure 4B2 (top 4 hits DSEA enrichment plots per cancer types)
plot_dsea_top_enrichment_plots_cancerspecific(
  dsea_results_df = dsea_cancer_specific_results,
  correlation_results_df = correlation_results_ANOVA,
  GDSC_age_res = GDSC_age_res,
  output_path = "Figures/Figure4B2.pdf",
  top_n = 2,
  minSize = 5,
  base_size = 26,
  title_size = 38,
  subtitle_size = 26,
  axis_title_size = 26,
  axis_text_size = 28,
  line_width = 2.4
)
