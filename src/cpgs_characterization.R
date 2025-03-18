##########################################################
###   Characterization of CpGs included in the model   ###
###   Author: Irene Fernández Rebollo                  ###
###   Date: 20/11/2024                                 ###
##########################################################

source("/Volumes/iorio/Irene/git_epiclock/src/utils.R")

# Load packages
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(missMethyl)
library(ggplot2)
library(VennDiagram)

# Load data
model_coefs <- read.csv("/Volumes/iorio/Irene/git_epiclock/res/model/model_coefs.csv")
model_coefs <- subset(model_coefs, s0 != 0) #4863  2
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_model_cpgs <- merge(ann, model_coefs, by.x = "row.names", by.y = "X")

# CpGs outside gene region
gene_region <- as.data.frame(ann_model_cpgs) %>% 
  mutate(group = ifelse(ann_model_cpgs$UCSC_RefGene_Name == "", "No gene", "Gene")) %>%
  group_by(group) %>%
  summarize(n = n()) %>%
  mutate(percentage = n / sum(n) * 100)

png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/cpgs_gene.png",
    width = 5, height = 5, units = 'in', res = 300)
ggplot(gene_region, aes(x = "", y = n, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c('deepskyblue3', 'grey80'),
                    labels = paste(gene_region$group, " (", round(gene_region$percentage, 1), "%)", sep = ""),
                    name = "") +
  theme_void(base_size = 20)
dev.off()

# Gene group representation
cpg_group <- as.data.frame(ann_model_cpgs) %>%
  mutate(group = strsplit(as.character(ann_model_cpgs$UCSC_RefGene_Group), ";")) %>% 
  unnest(group) %>%
  group_by(group) %>%
  summarize(n = n()) %>%
  mutate(percentage = n / sum(n) * 100)

gene_order <- c("TSS1500", "TSS200", "1stExon", "5'UTR", "Body", "3'UTR")
cpg_group$group <- factor(cpg_group$group, levels = gene_order)
cpg_group_sorted <- cpg_group %>%
  arrange(group)

png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/cpgs_region.png",
    width = 5, height = 5, units = 'in', res = 300)
ggplot(cpg_group_sorted, aes(x = "", y = n, fill = group)) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(length(unique(cpg_group_sorted$group)), "Set2"),
                    labels = paste(cpg_group_sorted$group, " (", round(cpg_group_sorted$percentage, 1), "%)", sep = ""),
                    name = "") +
  theme_void(base_size = 20)
dev.off()


# Enrichment analysis - No ranking
genes <- unique(unlist(strsplit(ann_model_cpgs[ann_model_cpgs$UCSC_RefGene_Name != "", "UCSC_RefGene_Name"], ";")))
writeLines(genes, "/Volumes/iorio/Irene/git_epiclock/res/resources/3337genes.txt")
enrichbp <- enrichGO(gene = genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_bp.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichbp, title = "GO Biological Process")
dev.off()
enrichmf <- enrichGO(gene = genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "MF")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_mf.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichmf, title = "GO Molecular Function")
dev.off()
enrichcc <- enrichGO(gene = genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "CC")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_cc.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichcc, title = "GO Cellular Component")
dev.
entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enrichkegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
                         organism = "hsa")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_kegg.png",
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
             filename = '/Volumes/iorio/Irene/git_epiclock/res/plots/cpgs_coefs_sign.png',
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
             filename = '/Volumes/iorio/Irene/git_epiclock/res/plots/genes_coefs_sign.png',
             output=TRUE)

enrichbp_mix <- enrichGO(gene = genes_mix, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_bp_mix.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichbp_mix, title = "GO Biological Process")
dev.off()
enrichmf_mix <- enrichGO(gene = genes_mix, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_mf_mix.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichmf_mix, title = "GO Molecular Function")
dev.off()
enrichcc_mix <- enrichGO(gene = genes_mix, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_cc_mix.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichcc_mix, title = "GO Cellular Component")
dev.off()
entrez_ids_mix <- bitr(genes_mix, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enrichkegg_mix <- enrichKEGG(gene = entrez_ids_mix$ENTREZID, organism = "hsa")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_kegg_mix.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichkegg_mix, title = "KEGG")
dev.off()

enrichbp_neg <- enrichGO(gene = genes_neg[!(genes_neg %in% genes_mix)], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_bp_neg.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichbp_neg, title = "GO Biological Process")
dev.off()
enrichmf_neg <- enrichGO(gene = genes_neg[!(genes_neg %in% genes_mix)], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_mf_neg.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichmf_neg, title = "GO Molecular Function")
dev.off()
enrichcc_neg <- enrichGO(gene = genes_neg[!(genes_neg %in% genes_mix)], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_cc_neg.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichcc_neg, title = "GO Cellular Component")
dev.off()
entrez_ids_neg <- bitr(genes_neg[!(genes_neg %in% genes_mix)], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enrichkegg_neg <- enrichKEGG(gene = entrez_ids_neg$ENTREZID, organism = "hsa")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_kegg_neg.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichkegg_neg, title = "KEGG")
dev.off()

enrichbp_pos <- enrichGO(gene = genes_pos[!(genes_pos %in% genes_mix)], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_bp_pos.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichbp_pos, title = "GO Biological Process")
dev.off()
enrichmf_pos <- enrichGO(gene = genes_pos[!(genes_pos %in% genes_mix)], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_mf_pos.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichmf_pos, title = "GO Molecular Function")
dev.off()
enrichcc_pos <- enrichGO(gene = genes_pos[!(genes_pos %in% genes_mix)], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_cc_pos.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichcc_pos, title = "GO Cellular Component")
dev.off()
entrez_ids_pos <- bitr(genes_pos[!(genes_pos %in% genes_mix)], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enrichkegg_pos <- enrichKEGG(gene = entrez_ids_pos$ENTREZID, organism = "hsa")
png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/enrichGO_kegg_pos.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichkegg_pos, title = "KEGG")
dev.off()


# Enrichment analysis - CpGs
go_cpgs_neg <- gometh(cpgs_neg, collection = "GO", array.type = "450K")
go_cpgs_pos <- gometh(cpgs_pos, collection = "GO", array.type = "450K")
min(go_cpgs_neg$FDR) # 1 - No significant results
min(go_cpgs_pos$FDR) # 1 - No significant results
kegg_cpgs_neg <- gometh(cpgs_neg, collection = "KEGG", array.type = "450K")
kegg_cpgs_pos <- gometh(cpgs_pos, collection = "KEGG", array.type = "450K")
min(kegg_cpgs_neg$FDR) # 1 - No significant results
min(kegg_cpgs_pos$FDR) # 0.2769831 - No significant results

