###############################################################
###   Impact of Missing Coefficients on Model Performance   ###
###   Author: Irene Fernández Rebollo                       ###
###   Date: 04/03/2025                                      ###
###############################################################

setwd("/group/iorio/Irene/epiclock_dev")
source("src/utils.R")
dir_results <- "results/TCGA/"
library(caret)

# Load data
cpgs_shared <- CpGshared() # CpGs shared between different arrays (450k, EPIC, EPICv2)
projects_tcga <- TCGA_projects() 
tcga_samples <- TCGA_samples()
tcga_bval <- TCGA_Bvalues()
# Create categorical age variable
tcga_samples$age_cat <- cut(tcga_samples$age_at_index, 
                            breaks = seq(0, max(tcga_samples$age_at_index, na.rm = TRUE) + 1, by = 10), 
                            right = FALSE, 
                            labels = paste0(seq(0, max(tcga_samples$age_at_index, na.rm = TRUE) - 1, by = 10),
                                            "-", seq(9, max(tcga_samples$age_at_index, na.rm = TRUE), by = 10)))
# Create variable combining project and age 
tcga_samples$project_age <- paste(tcga_samples$project, tcga_samples$age_cat, sep = "_")

# Divide data in model (tumor), no model (tumor), and normal
samples_normal <- tcga_samples[tcga_samples$sample_type == "Solid Tissue Normal",]
samples_tumor <- tcga_samples[tcga_samples$sample_type %in% c("Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood"), ]
samples_model <- samples_tumor[!(samples_tumor$patient %in% c(samples_normal$patient, samples_tumor[duplicated(samples_tumor$patient), "patient"])),]
samples_nomodel <- samples_tumor[!(samples_tumor$patient %in% samples_model$patient),]
bval_model <- tcga_bval[rownames(tcga_bval) %in% samples_model$barcode,]
bval_nomodel <- tcga_bval[rownames(tcga_bval) %in% samples_nomodel$barcode,]
bval_normal <- tcga_bval[rownames(tcga_bval) %in% samples_normal$barcode,]

# Train & Test split
set.seed(2025)
train_rows <- createDataPartition(samples_model[,"project_age"], p=0.8, list = FALSE)
train <- samples_model[train_rows, "barcode"]
test <- samples_model[-train_rows, "barcode"]

x.test <- bval_model[rownames(bval_model) %in% test, colnames(bval_model) %in% cpgs_shared]
y.test <- samples_model[samples_model$barcode %in% test, "age_at_index"]
x.normal <- bval_normal[, colnames(bval_normal) %in% cpgs_shared]
y.normal <- samples_normal[, "age_at_index"]
x.tumor <- bval_nomodel[, colnames(bval_nomodel) %in% cpgs_shared]
y.tumor <- samples_nomodel[, "age_at_index"]

# Load coefficients
coefs <- read.csv(paste0(dir_results, "model/model_coefs.csv"))
not0_coefs <- coefs[coefs[,2] != 0, 1]

# Performance by increasing the missing CpGs percentage
res_test <- c()
res_normal <- c()
res_tumor <- c()
missing_perc <- seq(0, 1, by = 0.05) # missing CpGs percentages
for(perc in missing_perc){
  for(i in 1:20){ # 20-CV
    missing_coefs <- sample(not0_coefs, ceiling(perc*length(not0_coefs)))
    pred_test <- AgePred_tcga(name = paste0("missing/test_", perc, "_", i, ".txt"),
                         coefs = coefs,
                         bval = t(x.test[, !(colnames(x.test) %in% missing_coefs)]))
    pred_normal <- AgePred_tcga(name = paste0("missing/normal_", perc, "_", i, ".txt"),
                           coefs = coefs,
                           bval = t(x.normal[, !(colnames(x.normal) %in% missing_coefs)]))
    pred_tumor <- AgePred_tcga(name = paste0("missing/tumor_", perc, "_", i, ".txt"),
                          coefs = coefs,
                          bval = t(x.tumor[, !(colnames(x.tumor) %in% missing_coefs)]))
    print("after AgePred")
    print(head(pred_test))
    r_test <- cor(y.test, pred_test[1,])
    r_normal <- cor(y.normal, pred_normal[1,])
    r_tumor <- cor(y.tumor, pred_tumor[1,])
    mae_test <- mean(abs(y.test - pred_test[1,]))
    mae_normal <- mean(abs(y.normal - pred_normal[1,]))
    mae_tumor <- mean(abs(y.tumor - pred_tumor[1,]))
    row_test <- c(perc, i, r_test, mae_test)
    row_normal <- c(perc, i, r_normal, mae_normal)
    row_tumor <- c(perc, i, r_tumor, mae_tumor)
    res_test <- rbind(res_test, row_test)
    res_normal <- rbind(res_normal, row_normal)
    res_tumor <- rbind(res_tumor, row_tumor)
  }
}
colnames(res_test) <- c("Perc_missing", "Iteration", "R", "MAE")
colnames(res_normal) <- c("Perc_missing", "Iteration", "R", "MAE")
colnames(res_tumor) <- c("Perc_missing", "Iteration", "R", "MAE")
write.csv(res_test, paste0(dir_results, "model/test_missing.csv"))
write.csv(res_normal, paste0(dir_results, "model/normal_missing.csv"))
write.csv(res_tumor, paste0(dir_results, "model/tumor_missing.csv"))

# Plot results
# # Load data
# res_test <- read.csv("/Volumes/iorio/Irene/epiclock_dev/results/TCGA/model/test_missing.csv")
# res_normal <- read.csv("/Volumes/iorio/Irene/epiclock_dev/results/TCGA/model/normal_missing.csv")
# res_tumor <- read.csv("/Volumes/iorio/Irene/epiclock_dev/results/TCGA/model/tumor_missing.csv")

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
png(filename = "/group/iorio/Irene/epiclock/plots/missing_coefficients.png",
    width = 12, height = 7, units = 'in', res = 600)
ggplot(res, aes(x = Perc_missing, y = mean,
                ymin = mean-sd, ymax = mean+sd, fill = data, color = data)) +
  geom_line() +
  geom_errorbar(width = 0.02) +
  geom_point(shape=21, color="black", size=3) +
  scale_fill_manual(values = color_type, name = "Data subset") +
  scale_color_manual(values = color_type, name = "Data subset") +
  theme_bw(base_size = 20) +
  xlab("Missing Coefficients (%)") + ylab("Pearson Correlation (R)") 
dev.off()
