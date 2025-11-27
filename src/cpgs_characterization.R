##########################################################
###   Characterization of CpGs included in the model   ###
###   Author: Irene Fernández Rebollo                  ###
###   Date: 20/11/2024                                 ###
##########################################################

setwd("/group/iorio/Irene/epiclock_dev")
source("src/utils.R")
dir_metadata <- "metadata/"

# Load packages
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(VennDiagram)
library(ggrepel)

# Load data
load(paste0(dir_metadata, "HorvathS2013.rda"))
horvath <- coefs
load(paste0(dir_metadata, "HannumG2013.rda"))
hannum <- coefs
load(paste0(dir_metadata, "HorvathS2018.rda"))
horvath2 <- coefs
load(paste0(dir_metadata, "LevineM2018.rda"))
levine <- coefs
model_coefs <- read.csv("results/TCGA/model/model_coefs.csv")
model_coefs <- subset(model_coefs, s0 != 0) #4863  2
cpgs_shared <- CpGshared()
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_horvath <- as.data.frame(merge(ann, horvath, by.x = "row.names", by.y = "Probe"))
ann_hannum <- as.data.frame(merge(ann, hannum, by.x = "row.names", by.y = "Probe"))
ann_horvath2 <- as.data.frame(merge(ann, horvath2, by.x = "row.names", by.y = "Probe"))
ann_levine <- as.data.frame(merge(ann, levine, by.x = "row.names", by.y = "Probe"))
ann_model_cpgs <- as.data.frame(merge(ann, model_coefs, by.x = "row.names", by.y = "X"))
ann_cpgs_shared <- as.data.frame(ann[rownames(ann) %in% cpgs_shared,])

# CpGs outside gene region
# Helper function to calculate gene region statistics
get_gene_region_stats <- function(data, source_name) {
  data %>%
    mutate(group = ifelse(UCSC_RefGene_Name == "", "No gene", "Gene")) %>%
    group_by(group) %>%
    summarize(n = n(), .groups = 'drop') %>%
    mutate(percentage = n / sum(n) * 100, source = source_name)
}

# Calculate statistics for both datasets
gene_region_horvath <- get_gene_region_stats(ann_horvath, "HorvathS2013")
gene_region_hannum <- get_gene_region_stats(ann_hannum, "HannumG2013")
gene_region_horvath2 <- get_gene_region_stats(ann_horvath2, "HorvathS2018")
gene_region_levine <- get_gene_region_stats(ann_levine, "LevineM2018")
gene_region <- get_gene_region_stats(ann_model_cpgs, "RebolloI2025")
gene_region_shared <- get_gene_region_stats(ann_cpgs_shared, "Illumina450k")

# Combine with Illumina array
combined_data <- bind_rows(gene_region,
                           gene_region_shared)

png(filename = "/Volumes/iorio/Irene/epiclock/plots/cpgs_gene.png",
    width = 10, height = 5, units = 'in', res = 600)
ggplot(combined_data, aes(x = source, y = percentage, fill = factor(group, levels = c("No gene", "Gene")))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('deepskyblue3', 'grey80'), name = "", breaks = c("Gene", "No gene")) +
  xlab("") + ylab("%CpGs") +
  coord_flip() +
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank())
dev.off()

# Combine with all clocks
combined_data <- bind_rows(gene_region_horvath,
                           gene_region_hannum,
                           gene_region_horvath2,
                           gene_region_levine,
                           gene_region)
combined_data$source <- factor(combined_data$source, levels = rev(c("HannumG2013", "HorvathS2013", "HorvathS2018", "LevineM2018", "RebolloI2025")))  # Replace with your desired order
png(filename = "/Volumes/iorio/Irene/epiclock/plots/cpgs_gene_clocks.png",
    width = 6, height = 4, units = 'in', res = 600)
ggplot(combined_data, aes(x = source, y = percentage, fill = factor(group, levels = c("Gene", "No gene")))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('deepskyblue3', 'grey80'), name = "", breaks = c("Gene", "No gene")) +
  scale_y_reverse(limits = c(100, 0)) +  # Reverse y-axis from 100 to 0
  xlab("Epigenetic Clocks") + ylab("%CpGs") +
  coord_flip() +
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank())
dev.off()

# png(filename = "results/TCGA/plots/cpgs_gene.png",
#     width = 5, height = 5, units = 'in', res = 300)
# ggplot(gene_region, aes(x = "", y = n, fill = group)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar(theta = "y") +
#   scale_fill_manual(values = c('deepskyblue3', 'grey80'),
#                     labels = paste(gene_region$group, " (", round(gene_region$percentage, 1), "%)", sep = ""),
#                     name = "") +
#   theme_void(base_size = 20)
# dev.off()

# Gene group representation
# Helper function to calculate gene region statistics
get_gene_group_stats <- function(data, source_name) {
  data %>%
    mutate(group = strsplit(as.character(UCSC_RefGene_Group), ";")) %>%
    unnest(group) %>%
    group_by(group) %>%
    summarize(n = n(), .groups = 'drop') %>%
    mutate(percentage = n / sum(n) * 100, source = source_name)
}

# Calculate statistics for both datasets
gene_group_horvath <- get_gene_group_stats(ann_horvath, "HorvathS2013")
gene_group_hannum <- get_gene_group_stats(ann_hannum, "HannumG2013")
gene_group_horvath2 <- get_gene_group_stats(ann_horvath2, "HorvathS2018")
gene_group_levine <- get_gene_group_stats(ann_levine, "LevineM2018")
gene_group <- get_gene_group_stats(ann_model_cpgs, "RebolloI2025")
gene_group_shared <- get_gene_group_stats(ann_cpgs_shared, "Illumina450k")

# Combine with Illumina array
combined_data <- bind_rows(gene_group,
                           gene_group_shared)

gene_order <- c("TSS1500", "TSS200", "1stExon", "5'UTR", "Body", "3'UTR")
combined_data$group <- factor(combined_data$group, levels = rev(gene_order))
combined_data_sorted <- combined_data %>%
  arrange(group)

png(filename = "/Volumes/iorio/Irene/epiclock/plots/cpgs_region.png",
    width = 6, height = 5, units = 'in', res = 600)
ggplot(combined_data, aes(x = source, y = percentage, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(length(gene_order), "Set2"),
                    breaks = gene_order,
                    name = "") +
  xlab("") + ylab("%CpGs") +
  coord_flip() +
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank())
dev.off()

# Combine with all clocks
combined_data <- bind_rows(gene_group_horvath,
                           gene_group_hannum,
                           gene_group_horvath2,
                           gene_group_levine,
                           gene_group)

gene_order <- c("TSS1500", "TSS200", "1stExon", "5'UTR", "Body", "3'UTR")
combined_data$group <- factor(combined_data$group, levels = rev(gene_order))
combined_data_sorted <- combined_data %>%
  arrange(group)

combined_data$source <- factor(combined_data$source, levels = rev(c("HannumG2013", "HorvathS2013", "HorvathS2018", "LevineM2018", "RebolloI2025")))  # Replace with your desired order
png(filename = "/Volumes/iorio/Irene/epiclock/plots/cpgs_region_clocks.png",
    width = 6, height = 4, units = 'in', res = 600)
ggplot(combined_data, aes(x = source, y = percentage, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(length(gene_order), "Set2"),
                    breaks = gene_order,
                    name = "") +
  xlab("Epigenetic Clocks") + ylab("%CpGs") +
  coord_flip() +
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank())
dev.off()
# png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/cpgs_region.png",
#     width = 5, height = 5, units = 'in', res = 300)
# ggplot(cpg_group_sorted, aes(x = "", y = n, fill = group)) +
#   geom_bar(stat = "identity") +
#   coord_polar(theta = "y") +
#   scale_fill_manual(values = RColorBrewer::brewer.pal(length(unique(cpg_group_sorted$group)), "Set2"),
#                     labels = paste(cpg_group_sorted$group, " (", round(cpg_group_sorted$percentage, 1), "%)", sep = ""),
#                     name = "") +
#   theme_void(base_size = 20)
# dev.off()

# Dotplot coefficients
gene_order <- c("TSS1500", "TSS200", "1stExon", "5'UTR", "Body", "3'UTR", "None")
model_coefs_ordered <- ann_model_cpgs %>%
  mutate(UCSC_RefGene_Group = ifelse(UCSC_RefGene_Group == "", "None", strsplit(UCSC_RefGene_Group, ";")),
         UCSC_RefGene_Name  = ifelse(UCSC_RefGene_Name == "", "None", strsplit(UCSC_RefGene_Name, ";"))) %>%
  unnest(c(UCSC_RefGene_Group, UCSC_RefGene_Name), keep_empty = TRUE) %>%
  distinct(Name, UCSC_RefGene_Name, UCSC_RefGene_Group, .keep_all = TRUE) %>%
  arrange(s0) %>%
  mutate(x = seq_along(s0),
         group = factor(UCSC_RefGene_Group, 
                        levels = rev(gene_order))) 
extreme_points <- model_coefs_ordered %>%
  filter(abs(s0) > 10) 

png(filename = "/Volumes/iorio/Irene/epiclock/plots/coefficients_region.png",
    width = 8, height = 5, units = 'in', res = 600)
ggplot(model_coefs_ordered, aes(x = x, y = s0, color = group)) +
  geom_jitter(width = 100, height = 0.2, size = 2) +
  scale_color_manual(values = RColorBrewer::brewer.pal(length(gene_order), "Set2"),
                     breaks = gene_order) +
  geom_text_repel(data = extreme_points, aes(x = x, y = s0, label = UCSC_RefGene_Name), 
                  size = 3, 
                  box.padding = 0.05,      # Increase padding around the labels
                  point.padding = 0.5,    # Increase space between the point and the label
                  max.overlaps = Inf,     # Remove overlap limit to give labels more flexibility
                  min.segment.length = 0, # Allow the labels to move freely
                  segment.size = 0.5, nudge_x = 400) + 
  ylab("Coefficient value") +
  theme_minimal(base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4)))
dev.off()


png(filename = "/Volumes/iorio/Irene/epiclock/plots/coefficients_region2.png",
    width = 7, height = 10, units = 'in', res = 600)
ggplot(model_coefs_ordered, aes(x = x, y = s0, color = group)) +
  geom_jitter(width = 500, height = 0.2, size = 2) +
  scale_color_manual(values = RColorBrewer::brewer.pal(length(gene_order), "Set2"),
                     breaks = gene_order) +
  geom_text_repel(data = extreme_points, aes(x = x, y = s0, label = UCSC_RefGene_Name), 
                  size = 4, 
                  box.padding = 0.2,      # Increase padding around the labels
                  point.padding = 0.9,    # Increase space between the point and the label
                  max.overlaps = Inf,     # Remove overlap limit to give labels more flexibility
                  min.segment.length = 0, # Allow the labels to move freely
                  segment.size = 0.5, nudge_x = 500) + 
  ylab("Coefficient value") +
  theme_minimal(base_size = 20) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4)))
dev.off()


# Enrichment analysis - No ranking
genes <- unique(unlist(strsplit(ann_model_cpgs[ann_model_cpgs$UCSC_RefGene_Name != "", "UCSC_RefGene_Name"], ";")))
writeLines(genes, "results/TCGA/resources/3337genes.txt")
enrichbp <- enrichGO(gene = genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP")
png(filename = "results/TCGA/plots/enrichGO_bp.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichbp, title = "GO Biological Process")
dev.off()
enrichmf <- enrichGO(gene = genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "MF")
png(filename = "results/TCGA/plots/enrichGO_mf.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichmf, title = "GO Molecular Function")
dev.off()
enrichcc <- enrichGO(gene = genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "CC")
png(filename = "results/TCGA/plots/enrichGO_cc.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichcc, title = "GO Cellular Component")
dev.
entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enrichkegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
                         organism = "hsa")
png(filename = "results/TCGA/plots/enrichGO_kegg.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichkegg, title = "KEGG")
dev.off()

# Enrichment analysis - coef sign
cpgs_neg <- model_coefs[model_coefs$s0 < 0, 1]
cpgs_pos <- model_coefs[model_coefs$s0 > 0, 1][-1]
venn.diagram(x = list(cpgs_neg, cpgs_pos),
             category.names = c("Negative" , "Positive"),
             col = c("red3", "green4"),
             fill = c(alpha("red3", 0.3), alpha("green4", 0.3)),
             cat.pos = c(-27, 27),
             cat.dist = c(0.055, 0.055),
             cat.fontfamily = "sans",
             fontfamily = "sans",
             imagetype = "png",
             height = 1500, 
             width = 1500,
             filename = 'results/TCGA/plots/cpgs_coefs_sign.png',
             output=TRUE)
ann_model_cpgs_neg <- ann_model_cpgs[ann_model_cpgs$Row.names %in% cpgs_neg,]
ann_model_cpgs_pos <- ann_model_cpgs[ann_model_cpgs$Row.names %in% cpgs_pos,]
genes_neg <- unique(unlist(strsplit(ann_model_cpgs_neg[ann_model_cpgs_neg$UCSC_RefGene_Name != "", "UCSC_RefGene_Name"], ";")))
genes_pos <- unique(unlist(strsplit(ann_model_cpgs_pos[ann_model_cpgs_pos$UCSC_RefGene_Name != "", "UCSC_RefGene_Name"], ";")))
genes_mix <- intersect(genes_neg, genes_pos)
venn.diagram(x = list(genes_neg, genes_pos),
             category.names = c("Negative" , "Positive"),
             col = c("red3", "green4"),
             fill = c(alpha("red3", 0.3), alpha("green4", 0.3)),
             cat.pos = c(-27, 27),
             cat.dist = c(0.055, 0.055),
             cat.fontfamily = "sans",
             fontfamily = "sans",
             imagetype = "png",
             height = 1500, 
             width = 1500,
             filename = 'results/TCGA/plots/genes_coefs_sign.png',
             output=TRUE)

enrichbp_mix <- enrichGO(gene = genes_mix, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
png(filename = "results/TCGA/plots/enrichGO_bp_mix.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichbp_mix, title = "GO Biological Process")
dev.off()
enrichmf_mix <- enrichGO(gene = genes_mix, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF")
png(filename = "results/TCGA/plots/enrichGO_mf_mix.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichmf_mix, title = "GO Molecular Function")
dev.off()
enrichcc_mix <- enrichGO(gene = genes_mix, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC")
png(filename = "results/TCGA/plots/enrichGO_cc_mix.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichcc_mix, title = "GO Cellular Component")
dev.off()
entrez_ids_mix <- bitr(genes_mix, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enrichkegg_mix <- enrichKEGG(gene = entrez_ids_mix$ENTREZID, organism = "hsa")
png(filename = "results/TCGA/plots/enrichGO_kegg_mix.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichkegg_mix, title = "KEGG")
dev.off()

enrichbp_neg <- enrichGO(gene = genes_neg[!(genes_neg %in% genes_mix)], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
png(filename = "results/TCGA/plots/enrichGO_bp_neg.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichbp_neg, title = "GO Biological Process")
dev.off()
enrichmf_neg <- enrichGO(gene = genes_neg[!(genes_neg %in% genes_mix)], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF")
png(filename = "results/TCGA/plots/enrichGO_mf_neg.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichmf_neg, title = "GO Molecular Function")
dev.off()
enrichcc_neg <- enrichGO(gene = genes_neg[!(genes_neg %in% genes_mix)], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC")
png(filename = "results/TCGA/plots/enrichGO_cc_neg.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichcc_neg, title = "GO Cellular Component")
dev.off()
entrez_ids_neg <- bitr(genes_neg[!(genes_neg %in% genes_mix)], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enrichkegg_neg <- enrichKEGG(gene = entrez_ids_neg$ENTREZID, organism = "hsa")
png(filename = "results/TCGA/plots/enrichGO_kegg_neg.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichkegg_neg, title = "KEGG")
dev.off()

enrichbp_pos <- enrichGO(gene = genes_pos[!(genes_pos %in% genes_mix)], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
png(filename = "results/TCGA/plots/enrichGO_bp_pos.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichbp_pos, title = "GO Biological Process")
dev.off()
enrichmf_pos <- enrichGO(gene = genes_pos[!(genes_pos %in% genes_mix)], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF")
png(filename = "results/TCGA/plots/enrichGO_mf_pos.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichmf_pos, title = "GO Molecular Function")
dev.off()
enrichcc_pos <- enrichGO(gene = genes_pos[!(genes_pos %in% genes_mix)], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC")
png(filename = "results/TCGA/plots/enrichGO_cc_pos.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichcc_pos, title = "GO Cellular Component")
dev.off()
entrez_ids_pos <- bitr(genes_pos[!(genes_pos %in% genes_mix)], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enrichkegg_pos <- enrichKEGG(gene = entrez_ids_pos$ENTREZID, organism = "hsa")
png(filename = "results/TCGA/plots/enrichGO_kegg_pos.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichkegg_pos, title = "KEGG")
dev.off()


# # Enrichment analysis - CpGs
# library(missMethyl)
# go_cpgs_neg <- gometh(cpgs_neg, collection = "GO", array.type = "450K")
# go_cpgs_pos <- gometh(cpgs_pos, collection = "GO", array.type = "450K")
# min(go_cpgs_neg$FDR) # 1 - No significant results
# min(go_cpgs_pos$FDR) # 1 - No significant results
# kegg_cpgs_neg <- gometh(cpgs_neg, collection = "KEGG", array.type = "450K")
# kegg_cpgs_pos <- gometh(cpgs_pos, collection = "KEGG", array.type = "450K")
# min(kegg_cpgs_neg$FDR) # 1 - No significant results
# min(kegg_cpgs_pos$FDR) # 0.2769831 - No significant results

