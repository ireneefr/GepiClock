#####################################################
###   Epigenetic Clocks Coefficients Comparison   ###
###   Author: Irene Fernández Rebollo             ###
###   Date: 06/06/2024                            ###
#####################################################

library(VennDiagram)
library(RColorBrewer)
library(scales)
library(UpSetR)
dir_metadata <- "/Volumes/iorio/Irene/git_epiclock/metadata/"

# Load data
load(paste0(dir_metadata, "HorvathS2013.rda"))
horvath <- coefs$Probe[-1]
load(paste0(dir_metadata, "HannumG2013.rda"))
hannum <- coefs$Probe
load(paste0(dir_metadata, "HorvathS2018.rda"))
horvath2 <- coefs$Probe[-1]
load(paste0(dir_metadata, "LevineM2018.rda"))
levine <- coefs$Probe[-1]
rebollo <- read.csv("/Volumes/iorio/Irene/git_epiclock/res/model/model_coefs.csv", row.names = 1)
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

# VennDiagram
myCol <- brewer.pal(5, "Pastel2")
venn.diagram(
  x = list(horvath,
           hannum,
           horvath,
           levine,
           rebollo),
  category.names = c("HorvathS2013",
                     "HannumG2013",
                     "HorvathS2018",
                     "LevineM2018",
                     "RebolloI2025"),
  filename = '/Volumes/iorio/Irene/git_epiclock/res/plots/venn_cpgs_comparison.png',
  fill = alpha(myCol, 0.7),
  col = myCol,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.dist = c(0.18,0.2,0.25,0.25,0.2),
  cat.pos = c(0,-20,-120,120,20),
  margin = 0.1,
  resolution = 300,
  output=TRUE
)
venn.diagram(
  x = list(gene_horvath,
           gene_hannum,
           gene_horvath,
           gene_levine,
           gene_rebollo),
  category.names = c("HorvathS2013",
                     "HannumG2013",
                     "HorvathS2018",
                     "LevineM2018",
                     "RebolloI2025"),
  filename = '/Volumes/iorio/Irene/git_epiclock/res/plots/venn_genes_comparison.png',
  fill = alpha(myCol, 0.7),
  col = myCol,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.dist = c(0.18,0.2,0.25,0.25,0.2),
  cat.pos = c(0,-20,-120,120,20),
  margin = 0.1,
  resolution = 300,
  output=TRUE
)

# Upset plot
cpg_lists <- list(HorvathS2013 = horvath,
                  HannumG2013 = hannum,
                  HorvathS2018 = horvath2,
                  LevineM2018 = levine,
                  RebolloI2025 = rebollo)
cpg_data_upset <- fromList(cpg_lists)

png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/upset_cpgs_comparison.png",
    width = 10, height = 6, units = 'in', res = 300)
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

gene_lists <- list(HorvathS2013 = gene_horvath,
                   HannumG2013 = gene_hannum,
                   HorvathS2018 = gene_horvath2,
                   LevineM2018 = gene_levine,
                   RebolloI2025 = gene_rebollo)
gene_data_upset <- fromList(gene_lists)

png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/upset_genes_comparison.png",
    width = 10, height = 6, units = 'in', res = 300)
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

