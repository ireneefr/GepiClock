#########################################################
###   Gene enrichment of CpGs included in the model   ###
###   Author: Irene Fernández Rebollo                 ###
###   Date: 22/01/2026                                ###
#########################################################

setwd("/Volumes/iorio/Irene/epiclock_dev")
source("src/utils.R")

# Load packages
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load data
model_coefs <- read.csv("results/TCGA/model/model_coefs.csv")
model_coefs <- subset(model_coefs, s0 != 0) #4863  2
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_model_cpgs <- as.data.frame(merge(ann, model_coefs, by.x = "row.names", by.y = "X"))

# Enrichment analysis - No ranking
genes <- unique(unlist(strsplit(ann_model_cpgs[ann_model_cpgs$UCSC_RefGene_Name != "", "UCSC_RefGene_Name"], ";")))
writeLines(genes, "results/TCGA/resources/3337genes.txt")
enrichbp <- enrichGO(gene = genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP")
enrichmf <- enrichGO(gene = genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "MF")
enrichcc <- enrichGO(gene = genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "CC")
entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enrichkegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
                         organism = "hsa")
png(filename = "Supplementary Figures/Supplementary_Figure_16A.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichbp, title = "GO Biological Process")
dev.off()
png(filename = "Supplementary Figures/Supplementary_Figure_16B.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichmf, title = "GO Molecular Function")
dev.off()
png(filename = "Supplementary Figures/Supplementary_Figure_16C.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichcc, title = "GO Cellular Component")
dev.off()
png(filename = "Supplementary Figures/Supplementary_Figure_16D.png",
    width = 7, height = 7, units = 'in', res = 300)
dotplot(enrichkegg, title = "KEGG")
dev.off()
