####################
###   Figure 3   ###
####################

setwd("/Volumes/iorio/Irene/epiclock_dev/")
source("src/utils.R")
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggbio)
library(tidyr)
library(dplyr)
library(philentropy) #for distance(X, method = "jensen-shannon")
library(pheatmap)

dir_metadata <- "metadata/"

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

# pdf("Figures/Figure3A.pdf", width = 8, height = 5)
# ggbio::autoplot(genomic_windows, aes(fill = CpG_model), layout = "karyogram") +
#   scale_fill_gradient(low = "white", high = "red", na.value = "white") +
#   labs(fill = 'CpGs/Mb') +
#   theme(panel.background = element_blank(),
#         strip.background = element_blank(),
#         axis.text = element_text(color = "black"),
#         text = element_text(size = 15),
#         legend.position = c(0.9, 0.5),
#         strip.placement.y = "outside",
#         strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5)) +
#   facet_grid(rows = vars(seqnames), switch = "y")
# dev.off()

pdf("Figures/Figure3A.pdf", width = 10, height = 7)
ggbio::autoplot(genomic_windows, aes(fill = CpG_model), layout = "karyogram") +
  scale_fill_gradient(low = "white", high = "red", na.value = "white",
                      guide = guide_colorbar(frame.colour = "black",
                                             frame.linewidth = 0.2,
                                             ticks.colour = "black")) +
  layout_karyogram(data = genomic_windows, aes(x = (start+end)/2, y = log_padj), geom = "line", color = "black", size = 0.2) +
  layout_karyogram(data = genomic_windows[genomic_windows$padj<=0.05,], aes(x = (start+end)/2, y = sig), geom = "point", shape = 7, color = "black", size = 0.25, stroke = 0.5) +
  labs(fill = 'CpGs/Mb') +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        panel.spacing.y = unit(0, "lines"))
dev.off()


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

# Heatmap for model CpG in chr (%)
pdf("Figures/Figure3A_1.pdf", width = 2, height = 8)
ggplot(df_chr, aes(x = "%CpGs in chr", y = chr, fill = perc_in_chr)) +
  geom_tile(color = "black", size = 0.2, width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", perc_in_chr)), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "slateblue1", name = "Model CpGs\nin chr (%)",
                      guide = guide_colorbar(
                        title.position = "bottom",   # title position: "top", "bottom", "left", "right"
                        title.hjust = 0.5, # center the title relative to bar
                        frame.colour = "black",
                        frame.linewidth = 0.2,
                        ticks.colour = "black")) +
  scale_y_discrete(limits = rev(levels(df_chr$chr))) +
  theme_void() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 14))
dev.off()

# Heatmap for chr CpG in model (%)
pdf("Figures/Figure3A_2.pdf", width = 2, height = 8)
ggplot(df_chr, aes(x = "%CpGs in model", y = chr, fill = perc_in_model)) +
  geom_tile(color = "black", size = 0.2, width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", perc_in_model)), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "aquamarine3", name = "Chr CpGs\nin model (%)",
                      guide = guide_colorbar(
                        title.position = "bottom",   # title position: "top", "bottom", "left", "right"
                        title.hjust = 0.5, # center the title relative to bar
                        frame.colour = "black",
                        frame.linewidth = 0.2,
                        ticks.colour = "black")) +
  scale_y_discrete(limits = rev(levels(df_chr$chr))) +
  theme_void() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 14))
dev.off()


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
gene_region_horvath <- get_gene_region_stats(ann_horvath, "Horvath\n2013")
gene_region_hannum <- get_gene_region_stats(ann_hannum, "Hannum")
gene_region_horvath2 <- get_gene_region_stats(ann_horvath2, "Horvath\n2018")
gene_region_levine <- get_gene_region_stats(ann_levine, "Levine")
gene_region <- get_gene_region_stats(ann_model_cpgs, "Rebollo")
gene_region_shared <- get_gene_region_stats(ann_cpgs_shared, "Illumina\narray")

# # Combine with Illumina array
# combined_data <- bind_rows(gene_region,
#                            gene_region_shared)
# 
# ggplot(combined_data, aes(x = source, y = percentage, fill = factor(group, levels = c("No gene", "Gene")))) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c('deepskyblue3', 'grey80'), name = "", breaks = c("Gene", "No gene")) +
#   xlab("") + ylab("%CpGs") +
#   coord_flip() +
#   theme_minimal(base_size = 20) +
#   theme(axis.text = element_text(color = "black"),
#         legend.position = "top",
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.x = element_blank())
# ggsave("Figures/Figure3B.pdf")

# Jensen-Shannon distance function
jsd <- function(X) {
  d <- distance(X, method = "jensen-shannon")
  m <- as.matrix(d)
  rownames(m) <- rownames(X)
  colnames(m) <- rownames(X)
  m
}

# Combine with all clocks
combined_data <- bind_rows(gene_region_horvath,
                           gene_region_hannum,
                           gene_region_horvath2,
                           gene_region_levine,
                           gene_region,
                           gene_region_shared)
combined_data$source <- factor(combined_data$source, 
                               levels = c("Illumina\narray", "Rebollo", "Levine",
                                          "Horvath\n2018", "Horvath\n2013", "Hannum"))
pdf("Figures/Figure3B.pdf", width = 8, height = 5)
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

df_wide <- combined_data %>%
  dplyr::select(source, group, percentage) %>%
  pivot_wider(names_from = group,
              values_from = percentage,
              values_fill = 0) %>%
  arrange(source)
X <- df_wide %>%
  dplyr::select(-source) %>%
  as.matrix() / 100
rownames(X) <- df_wide$source
jsd_mat <- jsd(X)
jsd_mat_lower <- jsd_mat

my_colors <- colorRampPalette(c("red3", "white"))(100)
pdf("Supplementary Figures/Supplementary_Figure_17A.pdf", width = 5, height = 5)
pheatmap(jsd_mat, 
         cluster_rows = FALSE,   # Removes row dendrograms
         cluster_cols = FALSE,   # Removes column dendrograms
         color = my_colors,      # Applies the simple 2-color gradient
         display_numbers = FALSE,
         angle_col = "45",
         border_color = "white",
         main = "Jensen-Shannon distance based on\nannotated gene regions")
dev.off()

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
gene_group_horvath <- get_gene_group_stats(ann_horvath, "Horvath\n2013")
gene_group_hannum <- get_gene_group_stats(ann_hannum, "Hannum")
gene_group_horvath2 <- get_gene_group_stats(ann_horvath2, "Horvath\n2018")
gene_group_levine <- get_gene_group_stats(ann_levine, "Levine")
gene_group <- get_gene_group_stats(ann_model_cpgs, "Rebollo")
gene_group_shared <- get_gene_group_stats(ann_cpgs_shared, "Illumina\narray")

# Combine with all clocks
combined_data <- bind_rows(gene_group_horvath,
                           gene_group_hannum,
                           gene_group_horvath2,
                           gene_group_levine,
                           gene_group,
                           gene_group_shared)
combined_data$source <- factor(combined_data$source, 
                               levels = c("Illumina\narray", "Rebollo", "Levine",
                                          "Horvath\n2018", "Horvath\n2013", "Hannum"))
gene_order <- c("TSS1500", "TSS200", "1stExon", "5'UTR", "Body", "3'UTR")
combined_data$group <- factor(combined_data$group, levels = rev(gene_order))
combined_data_sorted <- combined_data %>%
  arrange(group)
pdf("Figures/Figure2C.pdf", width = 8, height = 5.4)
ggplot(combined_data, aes(x = source, y = percentage, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(length(gene_order), "Set2"),
                    breaks = gene_order,
                    name = "CpGs within annotated\ngene regions") +
  xlab("") + ylab("%CpGs") +
  coord_flip() +
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank())
dev.off()

df_wide <- combined_data %>%
  dplyr::select(source, group, percentage) %>%
  pivot_wider(names_from = group,
              values_from = percentage,
              values_fill = 0) %>%
  arrange(source)
X <- df_wide %>%
  dplyr::select(-source) %>%
  as.matrix() / 100
rownames(X) <- df_wide$source
jsd_mat <- jsd(X)
jsd_mat_upper <- jsd_mat

my_colors <- colorRampPalette(c("red3", "white"))(100)
pdf("Supplementary Figures/Supplementary_Figure_17B.pdf", width = 5, height = 5)
pheatmap(jsd_mat,
         cluster_rows = FALSE,   # Removes row dendrograms
         cluster_cols = FALSE,   # Removes column dendrograms
         color = my_colors,      # Applies the simple 2-color gradient
         display_numbers = FALSE,
         angle_col = "45",
         border_color = "white",
         main = "Jensen-Shannon distance based on\nfunctional genomic annotation")
dev.off()
mean(jsd_mat[-2,2])













library(ComplexHeatmap)
library(circlize)
library(grid)

m <- jsd_mat
rownames(m) <- gsub("\n", " ", rownames(m))
colnames(m) <- gsub("\n", " ", colnames(m))
# two independent color mappings (adjust as you like)
col_lower <- colorRamp2(c(0, max(jsd_mat_lower, na.rm=TRUE)), c("red3", "black"))
col_upper <- colorRamp2(c(0, max(jsd_mat_upper, na.rm=TRUE)), c("blue3", "black"))

ht <- Heatmap(
  m,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = FALSE,
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_side = "top",
  row_names_side = "left",
  layer_fun = function(j, i, x, y, w, h, fill) {
    
    v <- m[cbind(i, j)]
    
    # ---- LOWER triangle ----
    selL <- i > j
    if (any(selL)) {
      fills <- col_lower(v[selL])
      grid.rect(x[selL], y[selL], w[selL], h[selL],
                gp = gpar(fill = fills, col = "white"))
      
      grid.text(sprintf("%.1e", v[selL]),
                x[selL], y[selL],
                gp = gpar(col = "white", fontsize = 10))
    }
    
    # ---- UPPER triangle ----
    selU <- i < j
    if (any(selU)) {
      fills <- col_upper(v[selU])
      grid.rect(x[selU], y[selU], w[selU], h[selU],
                gp = gpar(fill = fills, col = "white"))
      
      grid.text(sprintf("%.1e", v[selU]),
                x[selU], y[selU],
                gp = gpar(col = "white", fontsize = 10))
    }
    
    # ---- DIAGONAL (black) ----
    selD <- i == j
    if (any(selD)) {
      grid.rect(x[selD], y[selD], w[selD], h[selD],
                gp = gpar(fill = "white", col = "white"))
    }
  }
)


# two separate legends (one per triangle)
lgd_lower <- Legend(title = "Jensen-Shannon distance based on\nannotated gene regions",
                    col_fun = col_lower, direction = "horizontal", title_position = "topcenter")
lgd_upper <- Legend(title = "Jensen-Shannon distance based on\nfunctional genomic annotation",
                    col_fun = col_upper, direction = "horizontal", title_position = "topcenter")

pdf("Supplementary Figures/Supplementary_Figure_17A.pdf", width = 9, height = 7)
# Leave space: top = 12mm, bottom = 12mm
pushViewport(viewport(
  y = unit(0.5, "npc"),
  height = unit(1, "npc") - unit(40, "mm"),  # 12mm top + 12mm bottom
  just = "center"
))
draw(ht, show_heatmap_legend = FALSE, show_annotation_legend = FALSE, newpage = FALSE,
     padding = unit(c(4, 4, 4, 4), "mm"))
popViewport()
y_pos <- unit(1, "mm")

# shift distance from center
x_shift <- unit(4, "cm")

# Lower legend (left of center)
draw(lgd_lower,
     x = unit(0.55, "npc") - x_shift,
     y = y_pos,
     just = c("center", "bottom"))

# Upper legend (right of center)
draw(lgd_upper,
     x = unit(0.55, "npc") + x_shift,
     y = y_pos,
     just = c("center", "bottom"))
dev.off() 



# Plot distributions
vals_lower <- data.frame(jsd = jsd_mat_lower[upper.tri(jsd_mat_lower)])
vals_upper <- data.frame(jsd = jsd_mat_upper[upper.tri(jsd_mat_upper)])
pdf("Supplementary Figures/Supplementary_Figure_17B.pdf", width = 6, height = 4)
ggplot(vals_upper, aes(x = jsd)) +
  geom_density(color = "blue3", linewidth = 1, fill = "blue3", alpha = 0.3) +
  xlab("Jensen-Shannon distance based on\nfunctional genomic annotation") +
  theme_minimal(base_size = 15)
dev.off()
pdf("Supplementary Figures/Supplementary_Figure_17C.pdf", width = 6, height = 4)
ggplot(vals_lower, aes(x = jsd)) +
  geom_density(color = "red3", linewidth = 1, fill = "red3", alpha = 0.3) +
  xlab("Jensen-Shannon distance based on\nannotated gene regions") +
  theme_minimal(base_size = 15)
dev.off()
