##################################
###   Supplementary Figure 4   ###
##################################

source("src/utils.R")
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggbio)
library(tidyr)
library(dplyr)
dir_metadata <- "metadata/"

dir_metadata <- "metadata/"

# Load data
model_coefs <- read.csv("results/TCGA/model/model_coefs.csv")
model_coefs <- subset(model_coefs, s0 != 0) #4863  2
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_model_cpgs <- as.data.frame(merge(ann, model_coefs, by.x = "row.names", by.y = "X"))


# CpG models in karyogram
gr_model <- GRanges(seqnames = ann_model_cpgs$chr, 
                    ranges = IRanges(start = ann_model_cpgs$pos, end = ann_model_cpgs$pos),
                    coef = ann_model_cpgs$s0)
seqlevels(gr_model, pruning.mode = "coarse") <- paste0("chr", 1:22)
seqinfo(gr_model) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(gr_model)]

gr_all <- GRanges(seqnames = ann$chr, 
                  ranges = IRanges(start = ann$pos, end = ann$pos))
seqlevels(gr_all, pruning.mode = "coarse") <- paste0("chr", 1:22)
seqinfo(gr_all) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(gr_all)]

# CpGs per chromosome
df_chr <- data.frame(chr = paste0("chr", 1:22),
                     total_chr = sapply(paste0("chr", 1:22), function(x) length(gr_all[seqnames(gr_all) == x])),
                     model_chr = sapply(paste0("chr", 1:22), function(x) length(gr_model[seqnames(gr_model) == x])))

# Get percentages
total_all_cpgs <- sum(df_chr$total_chr)
total_model_cpgs <- sum(df_chr$model_chr)
df_chr <- df_chr %>%
  mutate(perc_in_chr = (model_chr / total_chr) * 100, # % model CpG in chr
         perc_in_model = (model_chr / total_model_cpgs) * 100, # % chr CpG in model
         perc_in_genome = (total_chr / total_all_cpgs) * 100) # % chr CpG in genome

# Ensure chromosomes are factors in the correct order
df_chr$chr <- factor(df_chr$chr, levels = paste0("chr", 1:22))

cpgs_shared <- CpGshared()
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_cpgs_shared <- as.data.frame(ann[rownames(ann) %in% cpgs_shared,])
df_chr <- as.data.frame(table(ann_cpgs_shared$chr))
df_chr <- df_chr[df_chr$Var1 %in% paste0("chr", 1:22),]
df_chr$Var1 <- factor(df_chr$Var1, levels = paste0("chr", 1:22))

# Calculate percentage of CpGs per chromosome
df_chr$perc_in_genome <- (df_chr$Freq / sum(df_chr$Freq)) * 100

# Histogram of CpGs across chromosomes in Illumina array
pdf("Supplementary Figures/Supplementary_Figure_4.pdf", width = 8, height = 5)
ggplot(df_chr, aes(x = Var1, y = perc_in_genome)) +
  geom_col(fill = "seagreen4") +
  theme_bw(base_size = 15) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(df_chr$Var1))) +
  xlab("") + ylab("Chromosomal CpG percentage") +
  theme(axis.text = element_text(color = "black"))
dev.off()
