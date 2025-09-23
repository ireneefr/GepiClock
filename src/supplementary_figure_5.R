##################################
###   Supplementary Figure 5   ###
##################################

library(UpSetR)
dir_metadata <- "metadata/"

# Load data
load(paste0(dir_metadata, "HorvathS2013.rda"))
horvath <- coefs$Probe[-1]
load(paste0(dir_metadata, "HannumG2013.rda"))
hannum <- coefs$Probe
load(paste0(dir_metadata, "HorvathS2018.rda"))
horvath2 <- coefs$Probe[-1]
load(paste0(dir_metadata, "LevineM2018.rda"))
levine <- coefs$Probe[-1]
rebollo <- read.csv("results/TCGA/model/model_coefs.csv", row.names = 1)
rebollo <- rownames(rebollo)[rebollo[,1] != 0][-1]

# Get genes
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_horvath <- ann[horvath,]
ann_hannum <- ann[hannum,]
ann_horvath2 <- ann[horvath2,]
ann_levine <- ann[levine,]
ann_rebollo <- ann[rebollo,]
gene_horvath <- unique(unlist(strsplit(ann_horvath$UCSC_RefGene_Name, ";")))
gene_hannum <- unique(unlist(strsplit(ann_hannum$UCSC_RefGene_Name, ";")))
gene_horvath2 <- unique(unlist(strsplit(ann_horvath2$UCSC_RefGene_Name, ";")))
gene_levine <- unique(unlist(strsplit(ann_levine$UCSC_RefGene_Name, ";")))
gene_rebollo <- unique(unlist(strsplit(ann_rebollo$UCSC_RefGene_Name, ";")))

# Upset plot
cpg_lists <- list('Horvath 2013' = horvath,
                  'Hannum' = hannum,
                  'Horvath 2018' = horvath2,
                  'Levine' = levine,
                  'Rebollo' = rebollo)
cpg_data_upset <- fromList(cpg_lists)

png("Supplementary Figures/Supplementary_Figure_5A.png",
    width = 10, height = 6, units = "in", res = 600)
upset(cpg_data_upset, order.by = "freq", nsets = 5,
      nintersects = 23,
      empty.intersections = "on",
      keep.order = TRUE,
      set_size.show = TRUE,
      mainbar.y.label = "Intersection of CpGs",
      sets.x.label = "CpG Coefficients", 
      set_size.scale_max = 6300,
      text.scale = c(1.75, 1.5, 1.5, 1.5, 1.75, 1.5))
dev.off()

gene_lists <- list('Horvath 2013' = gene_horvath,
                   'Hannum' = gene_hannum,
                   'Horvath 2018' = gene_horvath2,
                   'Levine' = gene_levine,
                   'Rebollo' = gene_rebollo)
gene_data_upset <- fromList(gene_lists)

png("Supplementary Figures/Supplementary_Figure_5B.png",
    width = 10, height = 6, units = "in", res = 600)
upset(gene_data_upset, order.by = "freq", nsets = 5,
      nintersects = 23,
      empty.intersections = "on",
      keep.order = TRUE,
      set_size.show = TRUE,
      mainbar.y.label = "Intersection of Genes",
      sets.x.label = "Genes", 
      set_size.scale_max = 6300,
      text.scale = c(1.75, 1.5, 1.5, 1.5, 1.75, 1.5))
dev.off()

