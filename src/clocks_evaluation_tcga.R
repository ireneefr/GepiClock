####################################################
###   Epigenetic Clocks Performance Comparison   ###
###   Author: Irene Fernández Rebollo            ###
###   Date: 25/02/2025                           ###
####################################################

setwd("/group/iorio/Irene/epiclock_dev")
source("src/utils.R")
library(caret)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
dir_results <- "results/TCGA/"

# Load data
tcga_samples <- TCGA_samples(dir_data = "/group/iorio/Irene/legacy/epiclock_old/data/")
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
project_pal <- read.csv("metadata/TCGAproject_palette.csv", col.names = c("project", "color"))


# Train & Test split
set.seed(2025)
train_rows <- createDataPartition(samples_model[,"project_age"], p=0.8, list = FALSE)
train <- samples_model[train_rows, "barcode"]
test <- samples_model[-train_rows, "barcode"]


# Combine information
tcga_samples_model <- tcga_samples_sub[tcga_samples_sub$barcode %in% c(train, test),]
tcga_samples_model[tcga_samples_model$barcode %in% train, "data_set"] <- "Primary Tumour (Train)"
tcga_samples_model[tcga_samples_model$barcode %in% test, "data_set"] <- "Primary Tumour (Test)"
tcga_samples_sub[tcga_samples_sub$barcode %in% test, "data_set"] <- "Primary Tumour (Test)"
tcga_samples_sub[tcga_samples_sub$barcode %in% samples_nomodel$barcode, "data_set"] <- "Primary Tumour (Others)"
tcga_samples_sub[tcga_samples_sub$barcode %in% samples_normal$barcode, "data_set"] <- "Normal Tissue"
tcga_samples_sub <- tcga_samples_sub[!is.na(tcga_samples_sub$data_set),]
tcga_samples_sub$sample_type <- factor(tcga_samples_sub$data_set, levels = c("Normal Tissue", "Primary Tumour (Test)", "Primary Tumour (Others)"))

tcga_samples_model <- merge(tcga_samples_model, pred_RebolloI2025, by.x = "barcode", by.y = "sample")
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

# Boxplot
# ggplot(tcga_samples_sub_plot, aes(x = clock, y = pred_age-age_at_index, fill = sample_type)) +
#   geom_boxplot() +
#   geom_hline(yintercept = 0, linetype = 2) +
#   scale_fill_manual(name = "Sample type", values = color_type) +
#   xlab("") + ylab("Predicted Age - Chronological Age") +
#   theme_minimal()

# Obtain correlations
correlation_by_group <- tcga_samples_sub_plot %>%
  group_by(Sample_type = sample_type, Clock = clock) %>%
  summarise(N = n(),
            R = cor.test(age_at_index, pred_age)$estimate,
            Pval = cor.test(age_at_index, pred_age)$p.value)

# Barplot
png(filename = "/Volumes/iorio/Irene/epiclock/plots/clocks_performance.png",
    width = 10, height = 7, units = 'in', res = 600)
ggplot(correlation_by_group, aes(x = Clock, y = R, fill = Sample_type)) +
  geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(name = "", values = color_type) +
  xlab("Epigenetic Clocks") + ylab("Pearson Correlation (R)") + ylim(c(0,1)) +
  scale_x_discrete(limits=rev) +
  theme_bw(base_size = 20) + 
  theme(axis.text = element_text(color = "black"),
        legend.position = "top") +
  coord_flip()
dev.off()

# Scatterplot
freq_sample_type <- as.data.frame(table(tcga_samples_sub$sample_type))
colnames(freq_sample_type) <- c("Sample_type", "Freq")
png(filename = "/Volumes/iorio/Irene/epiclock/plots/predvschron_RebolloI2025.png",
    width = 5, height = 10, units = 'in', res = 600)
tcga_samples_sub$Sample_type <- tcga_samples_sub$sample_type
ggplot(tcga_samples_sub, aes(x = age_at_index, y = RebolloI2025, fill = sample_type)) +
  geom_point(color = "black", shape = 21, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.75, aes(fill = NULL), show.legend = FALSE) +
  scale_fill_manual(name = "Data subset", values = color_type) +
  geom_text(data = freq_sample_type,
            aes(x = 5, y = 97, label = paste("N =", Freq)),  
            inherit.aes = FALSE, size = 5, hjust = 0) +
  geom_text(data = correlation_by_group[correlation_by_group$Clock == "RebolloI2025",],
            aes(x = 5, y = 90, label = paste("R =", round(R, 3))),  
            inherit.aes = FALSE, size = 5, hjust = 0) +
  geom_text(data = correlation_by_group[correlation_by_group$Clock == "RebolloI2025",],
            aes(x = 5, y = 83, label = paste("Pval =", formatC(Pval, format = "e", digits = 0))),  
            inherit.aes = FALSE, size = 5, hjust = 0) +
  facet_wrap(.~Sample_type, ncol = 1) +
  xlab("Chronological Age") + ylab("Predicted Age") +
  xlim(c(0, 100)) + ylim(c(0, 100)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color = "black"))
dev.off()

# # Obtain correlations by project - Not working :')
# correlation_by_group2 <- tcga_samples_sub_plot %>%
#   group_by(Sample_type = sample_type, Clock = clock, TCGA_project = project) %>%
#   summarise(N = n(),
#             R = ifelse(N >= 2, cor.test(age_at_index, pred_age)$estimate, NA),
#             Pval = ifelse(N >= 2, cor.test(age_at_index, pred_age)$p.value, NA))
# write.csv("/Volumes/iorio/Irene/git_epiclock/res/predvschron_RebolloI2025.png")
# 
# library(writexl)
# write_xlsx(list("Pearson correlations" = correlation_by_group,
#                 "Pearson correlations by project" = correlation_by_group2),
#            path = "/Volumes/iorio/Irene/git_epiclock/res/Pearson_correlations.xlsx")

# Correlation by project in test
project_correlation <- tcga_samples_model %>%
  group_by(project, data_set) %>%
  summarise(correlation = cor(age_at_index, RebolloI2025),
            N = n()) %>%
  arrange(correlation) %>%
  mutate(project = fct_reorder(project, correlation))
head(project_correlation)

p1 <- ggplot(project_correlation, aes(x = correlation, y = project, fill = data_set)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  theme_minimal(base_size = 15) +
  labs(x = "Correlation", y = NULL) +
  scale_fill_manual(values = c("steelblue1", "purple3")) +
  theme(legend.position = "top") +
  scale_x_reverse() +
  scale_y_discrete(position = "right")
p2 <- ggplot(project_correlation, aes(x = N, y = project, fill = data_set)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  theme_minimal(base_size = 15) +
  labs(x = "Sample Size (N)", y = NULL) +
  scale_fill_manual(values = c("steelblue1", "purple3")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
library(patchwork)
p1 + p2 + plot_layout(widths = c(1, 1))

ggplot(project_correlation, aes(x = N, y = correlation, color = project)) +
  geom_point()

library(dplyr)
library(ggplot2)
library(tidyr)

# Separate dataframes by group
N_data <- project_correlation %>% 
  select(project, data_set, N) %>% 
  rename(N_group = data_set, N_value = N)

corr_data <- project_correlation %>% 
  select(project, data_set, correlation) %>% 
  rename(corr_group = data_set, corr_value = correlation)

# Cross join
quad_data <- full_join(N_data, corr_data, by = "project")

# Create a label for each quadrant
quad_data <- quad_data %>%
  mutate(quadrant = paste0("N: ", N_group, "\nCorr: ", corr_group))
corr_labels <- quad_data %>%
  group_by(quadrant) %>%
  summarise(correlation = cor(N_value, corr_value, method = "pearson", use = "complete.obs")) %>%
  mutate(label = paste0("r = ", round(correlation, 2)))
ggplot(quad_data, aes(x = N_value, y = corr_value, color = project)) +
  geom_point(size = 3) +
  facet_wrap(~ quadrant, ncol = 2) +
  geom_text(data = corr_labels, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.2, inherit.aes = FALSE, size = 4) +
  labs(x = "N", y = "Correlation", title = "Quadrant Scatterplot by Group Combinations") +
  theme_minimal()

# Correlation by project in test
project_correlation <- tcga_samples_sub[tcga_samples_sub$sample_type == "Primary Tumour (Test)",] %>%
  group_by(project) %>%
  summarise(correlation = cor(age_at_index, RebolloI2025),
            N = n()) %>%
  arrange(correlation) %>%  # arrange by correlation
  mutate(project = factor(project, levels = project))

# Plot the correlation values in test
png(filename = "/Volumes/iorio/Irene/epiclock/plots/cor_test_project.png",
    width = 5, height = 10, units = 'in', res = 600)
ggplot(project_correlation, aes(y = project, x = correlation)) +
  geom_bar(stat = "identity") +  # stat = "identity" to use the precomputed correlation values
  theme_minimal(base_size = 15) +
  labs(y = "", x = "Pearson Correlation (R)")
dev.off()
# Plot the frequency values in test
png(filename = "/Volumes/iorio/Irene/epiclock/plots/freq_test_project.png",
    width = 5, height = 10, units = 'in', res = 600)
ggplot(project_correlation, aes(y = project, x = N)) +
  geom_bar(stat = "identity") +  # stat = "identity" to use the precomputed correlation values
  theme_minimal(base_size = 15) +
  labs(y = "", x = "Number of samples")
dev.off()

# Age distribution by project in test
tcga_samples_sub$project <- factor(tcga_samples_sub$project, 
                                   levels = levels(project_correlation$project))
png(filename = "/Volumes/iorio/Irene/epiclock/plots/age_test_project.png",
    width = 5, height = 10, units = 'in', res = 600)
ggplot(tcga_samples_sub[tcga_samples_sub$sample_type == "Primary Tumour (Test)",], aes(x = age_at_index, y = project)) +
  geom_boxplot() +
  theme_minimal(base_size = 15) +
  labs(y = "", x = "Chronological age")
dev.off()




# Correlation by project in test
project_correlation <- tcga_samples_sub[tcga_samples_sub$sample_type == "Normal Tissue",] %>%
  group_by(project) %>%
  summarise(correlation = cor(age_at_index, RebolloI2025),
            N = n()) %>%
  arrange(correlation) %>%  # arrange by correlation
  mutate(project = factor(project, levels = project))

# Plot the correlation values in test
png(filename = "/Volumes/iorio/Irene/epiclock/plots/cor_normal_project.png",
    width = 5, height = 10, units = 'in', res = 600)
ggplot(project_correlation, aes(y = project, x = correlation)) +
  geom_bar(stat = "identity") +  # stat = "identity" to use the precomputed correlation values
  theme_minimal(base_size = 15) +
  labs(y = "", x = "Pearson Correlation (R)")
dev.off()
# Plot the frequency values in test
png(filename = "/Volumes/iorio/Irene/epiclock/plots/freq_normal_project.png",
    width = 5, height = 10, units = 'in', res = 600)
ggplot(project_correlation, aes(y = project, x = N)) +
  geom_bar(stat = "identity") +  # stat = "identity" to use the precomputed correlation values
  theme_minimal(base_size = 15) +
  labs(y = "", x = "Number of samples")
dev.off()

# Age distribution by project in test
tcga_samples_sub$project <- factor(tcga_samples_sub$project, 
                                   levels = levels(project_correlation$project))
png(filename = "/Volumes/iorio/Irene/epiclock/plots/age_normal_project.png",
    width = 5, height = 10, units = 'in', res = 600)
ggplot(tcga_samples_sub[tcga_samples_sub$sample_type == "Normal Tissue",], aes(x = age_at_index, y = project)) +
  geom_boxplot() +
  theme_minimal(base_size = 15) +
  labs(y = "", x = "Chronological age")
dev.off()
