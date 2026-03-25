############################
###   FIGURE_5A         ###
###########################

# Generation of Figure 5A1 and 5A2

# ============================
library(tidyverse)
library(dplyr)
library(data.table)
library(fgsea)
library(writexl)
library(MutExMatSorting)
library(reshape2)
library(readxl)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# set working directory
setwd("/group/iorio/Alessandro.D/epiclock")

# load scripts
source("src/global_params.R")
source("src/utils.R")
source("src/plot_generation_helpers.R")

# load significant target enrichments (BP) in different cancer types
significant_enrichments_bp_cancer_spec <- read_csv(
  paste0(
    results_path,
    "cell_lines/cell_lines_gene_dependencies/gsea_per_cancer_type/01_significant_combined_BP_results.csv"
  ),
  show_col_types = FALSE
)

# check
stopifnot(all(c("ID", "Description", "NES", "p.adjust", "cancer_type") %in%
                colnames(significant_enrichments_bp_cancer_spec)))

# number of enrichments
for (cancer in names(table(significant_enrichments_bp_cancer_spec$cancer_type))) {
  print(paste0(
    "There are ",
    table(significant_enrichments_bp_cancer_spec$cancer_type)[cancer],
    " gene target enrichments (GOBP) in ",
    cancer
  ))
}

# load correlation results used to rank genes
correlation_results_ANOVA <- read_xlsx(
  paste0(
    results_path,
    "cell_lines/cell_lines_gene_dependencies/anova_cancer_specific_cumulative_correlation_results.xlsx"
  )
)

# check
print(head(correlation_results_ANOVA))
stopifnot(all(c("cancer_type", "Gene", "correlation") %in%
                colnames(correlation_results_ANOVA)))

# load GO Biological Process gene sets from MSigDB
msig_go_bp <- msigdbr(
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
)

# gene sets indexed by GO ID
go_bp_gene_sets_by_id <- split(msig_go_bp$gene_symbol, toupper(trimws(msig_go_bp$gs_exact_source)))
go_bp_gene_sets_by_id <- lapply(go_bp_gene_sets_by_id, unique)

# gene sets indexed by pathway description
go_bp_gene_sets_by_desc <- split(msig_go_bp$gene_symbol, tolower(trimws(msig_go_bp$gs_description)))
go_bp_gene_sets_by_desc <- lapply(go_bp_gene_sets_by_desc, unique)

# gene sets indexed by pathway name converted to readable text
msig_go_bp$gs_name_readable <- msig_go_bp$gs_name %>%
  gsub("^GOBP_", "", .) %>%
  gsub("_", " ", .) %>%
  tolower() %>%
  trimws()

go_bp_gene_sets_by_name <- split(msig_go_bp$gene_symbol, msig_go_bp$gs_name_readable)
go_bp_gene_sets_by_name <- lapply(go_bp_gene_sets_by_name, unique)

# select the top 3 pathways per cancer type
# this is the same logic used in the dotplot
significant_enrichments_bp_top3 <- significant_enrichments_bp_cancer_spec %>%
  filter(p.adjust < 0.05) %>%
  group_by(cancer_type) %>%
  arrange(desc(abs(NES)), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()

# plot figure 5A1 (dotplot for cancer specific GSEA, top 3 pathways per cancer type)
gsea_cancerspecific_dotplot(
  significant_enrichments_bp_cancer_spec,
  title = NULL,
  output_dir = "Figures/Figure5A1.pdf"
)

# plot figure 5A2 (GSEA curves for the same top 3 pathways per cancer type, absolute gene rank)
plot_gsea_top_enrichment_plots_cancerspecific(
  significant_results_bp_df = significant_enrichments_bp_cancer_spec,
  correlation_results_df = correlation_results_ANOVA,
  go_bp_gene_sets_by_id = go_bp_gene_sets_by_id,
  go_bp_gene_sets_by_desc = go_bp_gene_sets_by_desc,
  go_bp_gene_sets_by_name = go_bp_gene_sets_by_name,
  output_path = "Figures/gsea_curves/Figure5A2.pdf",
  selected_results_bp_df = significant_enrichments_bp_top3,
  pathways_to_plot = NULL,
  minSize = 5,
  p_adj_threshold = 0.05,
  rank_scale = "absolute",
  base_size = 26,
  title_size = 34,
  axis_title_size = 24,
  axis_text_size = 20,
  legend_title_size = 14,
  legend_text_size = 12,
  line_width = 2.4,
  tick_line_width = 0.9,
  box_line_width = 0.5
)
