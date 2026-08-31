##############################################
###   Mutational burden vs Age in TCGA     ###
###   Author: Irene Fernández Rebollo      ###
###   Date: 25/08/2026                     ###
##############################################

setwd("/group/iorio/Irene/epiclock_dev")
source("src/utils.R")
library(caret)
library(dplyr)
library(ggpubr)

# Load data
mut_load <- read_delim("metadata/mutational_burden/mutation-load_updated.txt")
head(mut_load)
cnv <- read.delim("metadata/cnv/seg_based_scores.tsv")
head(cnv)
tcga_samples <- TCGA_samples()

# Create categorical age variable
tcga_samples$age_cat <- cut(tcga_samples$age_at_index, 
                            breaks = seq(0, max(tcga_samples$age_at_index, na.rm = TRUE) + 1, by = 10), 
                            right = FALSE, 
                            labels = paste0(seq(0, max(tcga_samples$age_at_index, na.rm = TRUE) - 1, by = 10),
                                            "-", seq(9, max(tcga_samples$age_at_index, na.rm = TRUE), by = 10)))
# Create variable combining project and age 
tcga_samples$project_age <- paste(tcga_samples$project, tcga_samples$age_cat, sep = "_")
# Subset of data to plot
tcga_samples_sub <- tcga_samples[tcga_samples$sample_type %in% c("Solid Tissue Normal", "Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood"), ]

# Divide samples
samples_normal <- tcga_samples[tcga_samples$sample_type == "Solid Tissue Normal",]
samples_tumor <- tcga_samples[tcga_samples$sample_type %in% c("Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood"), ]
samples_model <- samples_tumor[!(samples_tumor$patient %in% c(samples_normal$patient, samples_tumor[duplicated(samples_tumor$patient), "patient"])),]
samples_nomodel <- samples_tumor[!(samples_tumor$patient %in% samples_model$patient),]

# Load predictions
pred_RebolloI2025 <-  read.csv("results/TCGA/predictions/pred_RebolloI2025.csv", col.names = c("sample", "RebolloI2025"))
pred_HorvathS2013 <- read.csv("results/TCGA/predictions/pred_HorvathS2013.csv", col.names = c("x", "sample", "HorvathS2013"))[,2:3]
pred_HannumG2013 <- read.csv("results/TCGA/predictions/pred_HannumG2013.csv", col.names = c("x", "sample", "HannumG2013"))[,2:3]
pred_HorvathS2018 <- read.csv("results/TCGA/predictions/pred_HorvathS2018.csv", col.names = c("x", "sample", "HorvathS2018"))[,2:3]
pred_LevineM2018 <- read.csv("results/TCGA/predictions/pred_LevineM2018.csv", col.names = c("x", "sample", "LevineM2018"))[,2:3]

# Train & Test split
set.seed(2025)
train_rows <- createDataPartition(samples_model[,"project_age"], p=0.8, list = FALSE)
train <- samples_model[train_rows, "barcode"]
test <- samples_model[-train_rows, "barcode"]

# Combine information
tcga_samples_sub[tcga_samples_sub$barcode %in% test, "data_set"] <- "Primary Tumour (Test)"
tcga_samples_sub[tcga_samples_sub$barcode %in% samples_nomodel$barcode, "data_set"] <- "Primary Tumour (Others)"
tcga_samples_sub[tcga_samples_sub$barcode %in% samples_normal$barcode, "data_set"] <- "Normal Tissue"
tcga_samples_sub <- tcga_samples_sub[!is.na(tcga_samples_sub$data_set),]
tcga_samples_sub$sample_type <- factor(tcga_samples_sub$data_set, levels = c("Primary Tumour (Test)", "Primary Tumour (Others)", "Normal Tissue"))

tcga_samples_sub <- merge(tcga_samples_sub, pred_RebolloI2025, by.x = "barcode", by.y = "sample")
tcga_samples_sub <- merge(tcga_samples_sub, pred_HorvathS2013, by.x = "barcode", by.y = "sample")
tcga_samples_sub <- merge(tcga_samples_sub, pred_HannumG2013, by.x = "barcode", by.y = "sample")
tcga_samples_sub <- merge(tcga_samples_sub, pred_HorvathS2018, by.x = "barcode", by.y = "sample")
tcga_samples_sub <- merge(tcga_samples_sub, pred_LevineM2018, by.x = "barcode", by.y = "sample")
tcga_samples_sub_plot <- tcga_samples_sub %>% 
  pivot_longer(c(HorvathS2013, HannumG2013, HorvathS2018, LevineM2018, RebolloI2025), names_to = "clock", values_to = "pred_age")
color_type <- c("Primary Tumour (Test)" = "steelblue1",
                "Primary Tumour (Others)" = "blueviolet",
                "Normal Tissue" = "chartreuse3")


tcga_samples_sub[1:5, 1:5]
mut_load[1:5,1:5]

# Adapt barcode to sample ID from mutatation load
tcga_samples_sub <- tcga_samples_sub %>%
  mutate(Tumor_Sample_ID = substr(barcode, 1, 15))
tcga_samples_sub_mut <- merge(tcga_samples_sub, mut_load, by = "Tumor_Sample_ID")
tcga_samples_sub_mut <- merge(tcga_samples_sub_mut, cnv, by.x = "Tumor_Sample_ID", by.y = "Sample")

# Plot mutational load vs age
pdf("results/mutational_burden/tcga_mut_load_age.pdf", height = 8, width = 12)
ggplot(tcga_samples_sub_mut, aes(x = age_at_index, y = `Non-silent per Mb`, colour = sample_type)) +
  geom_point() +
  stat_cor(aes(group = 1), method = "spearman", colour = "black",
           label.x.npc = "left", label.y.npc = "top", cor.coef.name = "rho") +
  scale_color_manual(values = color_type, name = "Sample type") +
  xlab("Chronological Age") +
  theme_minimal(base_size = 20)
dev.off()
pdf("results/mutational_burden/tcga_mut_load_gepiclock.pdf", height = 8, width = 12)
ggplot(tcga_samples_sub_mut, aes(x = RebolloI2025, y = `Non-silent per Mb`, colour = sample_type)) +
  geom_point() +
  stat_cor(aes(group = 1), method = "spearman", colour = "black",
           label.x.npc = "left", label.y.npc = "top", cor.coef.name = "rho") +
  scale_color_manual(values = color_type, name = "Sample type") +
  xlab("GepiClock Age") +
  theme_minimal(base_size = 20)
dev.off()

# Plot CNV burden vs age
pdf("results/mutational_burden/tcga_cnv_burden_age.pdf", height = 8, width = 12)
ggplot(tcga_samples_sub_mut, aes(x = age_at_index, y = frac_altered, colour = sample_type)) +
  geom_point() +
  stat_cor(aes(group = 1), method = "spearman", colour = "black",
           label.x.npc = "left", label.y.npc = "top", cor.coef.name = "rho") +
  scale_color_manual(values = color_type, name = "Sample type") +
  xlab("Chronological Age") +
  theme_minimal(base_size = 20)
dev.off()
pdf("results/mutational_burden/tcga_cnv_burden_gepiclock.pdf", height = 8, width = 12)
ggplot(tcga_samples_sub_mut, aes(x = RebolloI2025, y = frac_altered, colour = sample_type)) +
  geom_point() +
  stat_cor(aes(group = 1), method = "spearman", colour = "black",
           label.x.npc = "left", label.y.npc = "top", cor.coef.name = "rho") +
  scale_color_manual(values = color_type, name = "Sample type") +
  xlab("GepiClock Age") +
  theme_minimal(base_size = 20)
dev.off()
