############################
###   FIGURE_4A         ###
###########################

# Generation of Figure 4A (4A1,4A2 and 4A3) 

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
GDSC_age_res <- readRDS(paste0(results_path, "cell_lines/cell_lines_drug_sensitivity/cell_lines_per_drugs_pancancer.rds"))

# check
print(colnames(GDSC_age_res)[1:20])
stopifnot(all(c("DRUG_NAME", "PUTATIVE_TARGET") %in% colnames(GDSC_age_res)))

# load enrichment data
dsea_pan_cancer_results <- read_xlsx(paste0(results_path,
                                            "cell_lines/cell_lines_drug_sensitivity/dsea_pancancer/significant_enrichment_dsea_pancancer.xlsx"))

# load correlation results used to rank drugs
correlation_results_ANOVA <- read_xlsx(
  paste0(results_path, "cell_lines/cell_lines_drug_sensitivity/anova_pancancer_correlation_results.xlsx")
)
print(head(correlation_results_ANOVA))

# plot figure 4A1 (drug target enrichment)
ordered_targets <- dsea_pancancer_dotplot(dsea_result = dsea_pan_cancer_results, "Drug Target Enrichment by Predicted Group",
                                          output_dir = "Figures/Figure4A1.pdf")

# plot figure 4A2 (associated number of drugs)
plot_number_target_drug(
  significant_results_DSEA_df = dsea_pan_cancer_results,
  dsea_pan_cancer_results = dsea_pan_cancer_results,
  ordered_targets = ordered_targets,
  output_dir = "Figures/Figure4A2.pdf"
)

# plot figure 4A2 (associated drugs)
plot_target_drug_distribution(dsea_pan_cancer_results, ordered_targets, 
                              output_dir = "Figures/Figure4A2_alternative.pdf")

# plot figure 4A3 (top 4 hits DSEA enrichment plots)
plot_dsea_top_enrichment_plots(
  dsea_results_df = dsea_pan_cancer_results,
  correlation_results_df = correlation_results_ANOVA,
  GDSC_age_res = GDSC_age_res,
  output_path = "Figures/Figure4A3.pdf",
  top_n = 2,
  combine = FALSE,
  base_size = 26,
  title_size = 38,
  subtitle_size = 26,
  axis_title_size = 26,
  axis_text_size = 28,
  line_width = 2.4
)

