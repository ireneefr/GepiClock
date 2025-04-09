##########################################################
###   Genes represented by several CpGs in the model   ###
###   Author: Irene Fernández Rebollo                  ###
###   Date: 17/02/2025                                 ###
##########################################################

setwd("/group/iorio/Irene/epiclock_dev")
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)
library(ComplexHeatmap)
library(colorRamp2)
library(openxlsx)

# Load data
coefs <- read.csv("results/TCGA/model/model_coefs.csv", row.names = 1)
cpgs <- rownames(coefs)[coefs[,1] != 0][-1]

# Get genes
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_cpgs <- as.data.frame(ann[cpgs, c("chr", "Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group")])
ann_cpgs <- ann_cpgs[ann_cpgs$UCSC_RefGene_Name != "", ]
gene_cpgs <- unique(unlist(strsplit(ann_cpgs$UCSC_RefGene_Name, ";")))

# Split gene names and gene groups into separate rows
gene_long <- ann_cpgs %>%
  separate_rows(UCSC_RefGene_Name, UCSC_RefGene_Group, sep = ";") %>%
  distinct()

# Add coefficients info
gene_long <- merge(gene_long, coefs, by.x = "Name", by.y = "row.names")

# Calculate CpG count per gene separately
cpg_counts <- gene_long %>%
  distinct(Name, UCSC_RefGene_Name) %>%
  group_by(UCSC_RefGene_Name) %>%
  summarise(CpG_count = n(), .groups = "drop")

png(filename = "/Volumes/iorio/Irene/epiclock/plots/cpgs_by_gene.png",
    width = 5, height = 5, units = 'in', res = 600)
ggplot(cpg_counts, aes(x = as.factor(CpG_count))) +
  geom_bar() +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 4) +
  xlab("CpGs/Gene") +
  ylab("Number of Genes") +
  theme_minimal(base_size = 20)
dev.off()

# Merge the counts back into the original dataframe without losing rows
gene_long <- gene_long %>%
  left_join(cpg_counts, by = "UCSC_RefGene_Name")

# Get ranges for lollipop plot
gene_ranges <- gene_long %>%
  group_by(UCSC_RefGene_Name) %>%
  summarize(min_coeff = min(s0),
            max_coeff = max(s0),
            CpG_count = unique(CpG_count))

# Plot
png(filename = "/Volumes/iorio/Irene/epiclock/plots/overrepresented_genes.png",
    width = 10, height = 6, units = 'in', res = 600)
ggplot(gene_long[gene_long$CpG_count > 5,], aes(x = s0, y = reorder(UCSC_RefGene_Name, CpG_count))) +
  facet_wrap(~ CpG_count, scales = "free_y", labeller = labeller(CpG_count = function(x) paste(x, "CpGs/Gene"))) +
  geom_segment(data = gene_ranges[gene_ranges$CpG_count > 5,],
               aes(x = min_coeff, xend = max_coeff, y = UCSC_RefGene_Name, yend = UCSC_RefGene_Name),
               color = "black", size = 0.75, alpha = 0.6) +
  geom_point(aes(color = UCSC_RefGene_Group), size = 3, alpha = 0.75,
             position = position_jitter(width = 0.1, height = 0)) +
  labs(x = "Coefficient",
       y = "Genes") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom")
dev.off()

# Heatmap
gene_long$Name_Group <- paste0(gene_long$UCSC_RefGene_Group, gene_long$Name)
# Count CpGs per gene and arrange
gene_order <- gene_long %>%
  arrange(desc(CpG_count)) %>%
  pull(UCSC_RefGene_Name) %>%
  unique()

# Order CpGs by RefGene group
order <- c("TSS1500", "TSS200", "1stExon", "5'UTR", "Body", "3'UTR")
cpg_order <- gene_long %>%
  mutate(UCSC_RefGene_Group = factor(UCSC_RefGene_Group, levels = order)) %>%
  arrange(UCSC_RefGene_Group) %>%
  pull(Name_Group)  %>%
  unique()

# Pivot to wide matrix: genes as rows, CpGs as columns
heatmap_matrix <- gene_long[gene_long$CpG_count > 4,] %>%
  select(UCSC_RefGene_Name, Name_Group, s0) %>%
  pivot_wider(names_from = Name_Group, values_from = s0) %>%
  column_to_rownames("UCSC_RefGene_Name")

# Reorder rows and columns
heatmap_matrix <- heatmap_matrix[gene_order[gene_order %in% rownames(heatmap_matrix)], 
                                 cpg_order[cpg_order %in% colnames(heatmap_matrix)]]

col_fun <- colorRamp2(c(min(heatmap_matrix, na.rm = TRUE),
                        max(heatmap_matrix, na.rm = TRUE)),
                      c("blue", "red"))
Heatmap(heatmap_matrix,
        name = "s0",
        col = col_fun,
        na_col = "white",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        column_names_rot = 90,
        heatmap_legend_param = list(title = "s0"))
# Save info
write.xlsx(list("Gene_counts" = cpg_counts[order(cpg_counts$CpG_count),],
                "Gene_counts_annotation" = gene_long),
                "results/TCGA/overrepresented_genes.xlsx",
                colNames = TRUE, rowNames = FALSE)
