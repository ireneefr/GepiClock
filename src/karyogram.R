###################################################
###   Karyogram of CpGs included in the model   ###
###   Author: Irene Fernández Rebollo           ###
###   Date: 20/11/2024                          ###
###################################################

setwd("/group/iorio/Irene/epiclock_dev")
source("src/utils.R")

# Load packages
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggbio)
library(tidyr)
library(dplyr)

# Load data
model_coefs <- read.csv("results/TCGA/model/model_coefs.csv")
model_coefs <- subset(model_coefs, s0 != 0) #4863  2
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_model_cpgs <- merge(ann, model_coefs, by.x = "row.names", by.y = "X")

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

# Create windows 1Mb
window_size <- 1e6
genomic_windows <- tileGenome(seqinfo(gr_model), tilewidth = window_size, cut.last.tile.in.chrom = TRUE) #try gr_all

# Mean coefficient by window
overlaps <- findOverlaps(genomic_windows, gr_model, select = "all")
mean_coef_per_window <- sapply(split(gr_model$coef[subjectHits(overlaps)], queryHits(overlaps)), mean, na.rm = TRUE)
genomic_windows$mean_coef <- NA  # Initialize with NA
genomic_windows$mean_coef[as.numeric(names(mean_coef_per_window))] <- mean_coef_per_window

# Compute hypergeometric test by chromosome
results_list <- list()
for(chr in unique(seqnames(genomic_windows))){
  genomic_windows_chr <- genomic_windows[seqnames(genomic_windows) == chr,]
  genomic_windows_chr$CpG_model <- countOverlaps(genomic_windows_chr, gr_model)
  genomic_windows_chr$CpG_chr <- countOverlaps(genomic_windows_chr, gr_all)
  
  # Total number of CpGs in the 450K array
  M <- length(gr_all[seqnames(gr_all) == chr,])
  # Total number of CpGs in the model
  N <- length(gr_model[seqnames(gr_model) == chr,])
  
  pval <- mapply(function(x, K) {
    phyper(q = x - 1,  # successes in sample minus 1 (P(X >= x))
           m = N,      # total successes in population (model CpGs)
           n = M - N,  # total failures in population (non-model CpGs)
           k = K,      # sample size (CpGs in window)
           lower.tail = FALSE) # P-value for enrichment
  }, genomic_windows_chr$CpG_model, genomic_windows_chr$CpG_chr)
  padj <- p.adjust(pval, method = "BH")
  log_padj <- log2(padj)
  
  # Store results as a data frame
  results_list[[chr]] <- data.frame(CpG_model = genomic_windows_chr$CpG_model,
                                    CpG_chr = genomic_windows_chr$CpG_chr,
                                    pval = pval,
                                    padj = padj,
                                    log_padj = log_padj)
}
df_test <- do.call(rbind, results_list)
mcols(genomic_windows) <- cbind(mcols(genomic_windows), df_test)

# Plot results
genomic_windows$sig <- ifelse(genomic_windows$padj<=0.05, 10, NA)
genomic_windows$ypos <- 1

# png(filename = "/Volumes/iorio/Irene/git_epiclock/res/plots/cpgs_karyogram.png",
#     width = 10, height = 6, units = 'in', res = 300)
# ggbio::autoplot(genomic_windows, aes(fill = CpG_model), layout = "karyogram") +
#   scale_fill_gradient(low = "white", high = "red", na.value = "white") +
#   labs(fill = 'CpGs/Mb') +
#   theme(panel.background = element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 10, face = "bold"))
# dev.off()

png(filename = "results/TCGA/plots/cpgs_karyogram_hypergeometric.png",
    width = 10, height = 6, units = 'in', res = 300)
ggbio::autoplot(genomic_windows, aes(fill = CpG_model), layout = "karyogram") +
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  layout_karyogram(data = genomic_windows, aes(x = (start+end)/2, y = log_padj), geom = "line", color = "black", size = 0.4) +
  layout_karyogram(data = genomic_windows[genomic_windows$padj<=0.05,], aes(x = (start+end)/2, y = sig), geom = "point", shape = 8, color = "black", size = 2, stroke = 1) +
  labs(fill = 'CpGs/Mb') +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        panel.spacing.y = unit(0, "lines"))
dev.off()

png(filename = "results/TCGA/plots/cpgs_karyogram_hypergeometric_coefs.png",
    width = 10, height = 6, units = 'in', res = 300)
ggbio::autoplot(genomic_windows, aes(fill = mean_coef), layout = "karyogram") +
  scale_fill_gradient2(low = "orange", mid = "white", high = "purple", na.value = "white") +
  layout_karyogram(data = genomic_windows, aes(x = (start+end)/2, y = log_padj), geom = "line", color = "black", size = 0.4) +
  layout_karyogram(data = genomic_windows[genomic_windows$padj<=0.05,], aes(x = (start+end)/2, y = sig), geom = "point", shape = 8, color = "black", size = 2, stroke = 1) +
  labs(fill = 'Mean coefficient\nCpG/Mb') +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        panel.spacing.y = unit(0, "lines"))
dev.off()


