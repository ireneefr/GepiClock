####################
###   Figure 2   ###
####################

source("src/utils.R")
library(caret)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Load data
tcga_samples <- TCGA_samples(dir_data = "/Volumes/iorio/Irene/legacy/epiclock_old/data/")
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
tcga_samples_sub$sample_type <- factor(tcga_samples_sub$data_set, levels = c("Normal Tissue", "Primary Tumour (Test)", "Primary Tumour (Others)"))

tcga_samples_sub <- merge(tcga_samples_sub, pred_RebolloI2025, by.x = "barcode", by.y = "sample")
tcga_samples_sub <- merge(tcga_samples_sub, pred_HorvathS2013, by.x = "barcode", by.y = "sample")
tcga_samples_sub <- merge(tcga_samples_sub, pred_HannumG2013, by.x = "barcode", by.y = "sample")
tcga_samples_sub <- merge(tcga_samples_sub, pred_HorvathS2018, by.x = "barcode", by.y = "sample")
tcga_samples_sub <- merge(tcga_samples_sub, pred_LevineM2018, by.x = "barcode", by.y = "sample")
tcga_samples_sub_plot <- tcga_samples_sub %>% 
  pivot_longer(c(HorvathS2013, HannumG2013, HorvathS2018, LevineM2018, RebolloI2025), names_to = "clock", values_to = "pred_age")
color_type <- c("Primary Tumour (Test)" = "steelblue1",
                "Primary Tumour (Others)" = "steelblue3",
                "Normal Tissue" = "chartreuse3")

# Obtain correlations
correlation_by_group <- tcga_samples_sub_plot %>%
  group_by(Sample_type = sample_type, Clock = clock) %>%
  summarise(N = n(),
            R = cor.test(age_at_index, pred_age)$estimate,
            Pval = cor.test(age_at_index, pred_age)$p.value)

# Barplot
ggplot(correlation_by_group, aes(x = Clock, y = R, fill = Sample_type)) +
  geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(name = "", values = color_type,
                    labels = c("Normal", "Tumour (test)", "Tumour (other)")) +
  xlab("") + ylab("Pearson Correlation (R)") + ylim(c(0,1)) +
  scale_x_discrete(limits=rev, labels = c("HannumG2013" = "Hannum (2013)",
                                          "HorvathS2013" = "Horvath (2013)",
                                          "HorvathS2018" = "Horvath (2018)",
                                          "LevineM2018" = "Levine (2018)",
                                          "RebolloI2025" = "Rebollo")) +
  theme_bw(base_size = 20) + 
  theme(axis.text = element_text(color = "black"),
        legend.position = "bottom") +
  coord_flip()
ggsave("Figures/Figure2B.pdf", dpi = 600, width = 8, height = 6)

# Scatterplot
freq_sample_type <- as.data.frame(table(tcga_samples_sub$sample_type))
colnames(freq_sample_type) <- c("Sample_type", "Freq")
tcga_samples_sub$Sample_type <- tcga_samples_sub$sample_type
ggplot(tcga_samples_sub, aes(x = age_at_index, y = RebolloI2025, fill = sample_type)) +
  geom_point(color = "black", shape = 21, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.75, aes(fill = NULL), show.legend = FALSE) +
  scale_fill_manual(name = "Data subset", values = color_type) +
  geom_text(data = freq_sample_type,
            aes(x = 5, y = 97, label = paste("N =", Freq)),  
            inherit.aes = FALSE, size = 4, hjust = 0) +
  geom_text(data = correlation_by_group[correlation_by_group$Clock == "RebolloI2025",],
            aes(x = 5, y = 90, label = paste("R =", round(R, 3))),  
            inherit.aes = FALSE, size = 4, hjust = 0) +
  geom_text(data = correlation_by_group[correlation_by_group$Clock == "RebolloI2025",],
            aes(x = 5, y = 83, label = paste("Pval =", formatC(Pval, format = "e", digits = 0))),  
            inherit.aes = FALSE, size = 4, hjust = 0) +
  facet_wrap(.~Sample_type, ncol = 1, labeller = as_labeller(c("Normal Tissue" = "Normal",
                                                               "Primary Tumour (Test)" = "Tumour (test)",
                                                               "Primary Tumour (Others)" = "Tumour (other)"))) +
  xlab("Chronological Age") + ylab("Predicted Age") +
  xlim(c(0, 100)) + ylim(c(0, 100)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        axis.text = element_text(color = "black"))
ggsave("Figures/Figure2A.pdf", dpi = 600, width = 5, height = 10)


# Missing coefficients plot
# Load data
res_test <- read.csv("/Volumes/iorio/Irene/epiclock_dev/results/TCGA/model/test_missing.csv")
res_normal <- read.csv("/Volumes/iorio/Irene/epiclock_dev/results/TCGA/model/normal_missing.csv")
res_tumor <- read.csv("/Volumes/iorio/Irene/epiclock_dev/results/TCGA/model/tumor_missing.csv")

# Summarize information
res_test_summary <- as.data.frame(res_test) %>%
  group_by(Perc_missing) %>%
  summarise(
    sd = as.numeric(sd(R, na.rm = TRUE)),
    mean = as.numeric(mean(R)),
    mean_mae = as.numeric(mean(MAE)),
    data = "Primary Tumour (Test)"
  )
res_normal_summary <- as.data.frame(res_normal) %>%
  group_by(Perc_missing) %>%
  summarise(
    sd = as.numeric(sd(R, na.rm = TRUE)),
    mean = as.numeric(mean(R)),
    mean_mae = as.numeric(mean(MAE)),
    data = "Normal Tissue"
  )
res_tumor_summary <- as.data.frame(res_tumor) %>%
  group_by(Perc_missing) %>%
  summarise(
    sd = as.numeric(sd(R, na.rm = TRUE)),
    mean = as.numeric(mean(R)),
    mean_mae = as.numeric(mean(MAE)),
    data = "Primary Tumour (Others)"
  )
res <- rbind(res_test_summary, res_normal_summary, res_tumor_summary)

# Plot
color_type <- c("Primary Tumour (Test)" = "steelblue1",
                "Primary Tumour (Others)" = "steelblue3",
                "Normal Tissue" = "chartreuse3")
ggplot(res, aes(x = Perc_missing, y = mean,
                ymin = mean-sd, ymax = mean+sd, fill = data, color = data)) +
  geom_line(linewidth = 1) +
  geom_errorbar(width = 0.02) +
  geom_point(shape=21, color="black", size=3) +
  scale_fill_manual(values = color_type, name = "Data subset") +
  scale_color_manual(values = color_type, name = "Data subset") +
  xlab("Missing Coefficients (%)") + ylab("Pearson Correlation (R)") +
  theme_bw(base_size = 20) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"))
ggsave("Figures/Figure2C.pdf", dpi = 600, width = 8, height = 6)
