##################################
###   Supplementary Figure 2   ###
##################################

source("src/utils.R")
library(VennDiagram)
dir_arrays <- "metadata/Illumina_arrays/"

# Load Illumina methylation arrays
ann_epicv2 <- data.table::fread(paste0(dir_arrays, "EPIC-8v2-0_A1.csv"), header = TRUE, fill = TRUE)[c(-1,-2,-3,-4,-5,-6,-7),]
ann_epic <- data.table::fread(paste0(dir_arrays, "infinium-methylationepic-v-1-0-b5-manifest-file.csv"), header = TRUE, fill = TRUE)[c(-1,-2,-3,-4,-5,-6,-7),]
ann_450k <-  data.table::fread(paste0(dir_arrays, "humanmethylation450_15017482_v1-2.csv"), header = TRUE, fill = TRUE)[c(-1,-2,-3,-4,-5,-6,-7),]

# Remove control CpGs
cpgs_epicv2 <- ann_epicv2$Inc.[1:which(ann_epicv2$Illumina == "[Controls]")-1]
cpgs_epic <- ann_epic$Inc.[1:which(ann_epic$Illumina == "[Controls]")-1]
cpgs_450k <- ann_450k$Inc.[1:which(ann_450k$Illumina == "[Controls]")-1]

# Get Supplementary Figure 2
venn.diagram(list(`Infinium\nHumanMethylation450\nBeadChip` = cpgs_450k,
                  `Infinium\nMethylationEPIC\nBeadChip` = cpgs_epic,
                  `Infinium\nMethylationEPIC v2.0\nBeadChip` = cpgs_epicv2),
             filename = "Supplementary Figures/Supplementary_Figure_2.pdf",
             fill = c("#1f77b4", "#ff7f0e", "#2ca02c"),
             fontfamily = "sans",
             cat.fontfamily = "sans",
             cat.pos = c(-35, 35, 180),
             cat.dist = c(0.13, 0.13, 0.13),
             margin = 0.05, 
             fontface = c("plain", "plain", "plain", "plain", "bold", "plain", "plain"))

