################################################################
###   Impact of Reducing Coefficients on Model Performance   ###
###   Author: Irene Fernández Rebollo                        ###
###   Date: 15/04/2025                                       ###
################################################################

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
coefs <- coefs[-1,]
not0_coefs <- coefs[coefs[,2] != 0, ]

# Order CpGs by the coefficient value
not0_coefs_dec <- not0_coefs[order(abs(not0_coefs[,2]), decreasing = TRUE), ]
not0_coefs_inc <- not0_coefs[order(abs(not0_coefs[,2]), decreasing = FALSE), ]

# # Performance by reducing the CpGs percentage
# res_test_dec <- c()
# res_normal_dec <- c()
# res_tumor_dec <- c()
# res_test_inc <- c()
# res_normal_inc <- c()
# res_tumor_inc <- c()
# reducing_perc <- seq(1, 0, by = -0.05) # reducing CpGs percentages
# for(perc in reducing_perc){
#   num_to_keep <- ceiling(perc*nrow(not0_coefs))
#   # Coefficients to keep (the most extreme)
#   reducing_coefs_dec <- not0_coefs_dec[seq_len(num_to_keep), 1]
#   # Coefficients to keep (the least extreme)
#   reducing_coefs_inc <- not0_coefs_inc[seq_len(num_to_keep), 1]
#   pred_test_dec <- AgePred_tcga(name = paste0("reducing/test_dec_", perc, ".txt"),
#                                 coefs = coefs,
#                                 bval = t(x.test[, !(colnames(x.test) %in% reducing_coefs_dec)]))
#   pred_normal_dec <- AgePred_tcga(name = paste0("reducing/normal_dec_", perc, ".txt"),
#                                   coefs = coefs,
#                                   bval = t(x.normal[, !(colnames(x.normal) %in% reducing_coefs_dec)]))
#   pred_tumor_dec <- AgePred_tcga(name = paste0("reducing/tumor_dec_", perc, ".txt"),
#                                  coefs = coefs,
#                                  bval = t(x.tumor[, !(colnames(x.tumor) %in% reducing_coefs_dec)]))
#   pred_test_inc <- AgePred_tcga(name = paste0("reducing/test_inc_", perc, ".txt"),
#                                 coefs = coefs,
#                                 bval = t(x.test[, !(colnames(x.test) %in% reducing_coefs_inc)]))
#   pred_normal_inc <- AgePred_tcga(name = paste0("reducing/normal_inc_", perc, ".txt"),
#                                   coefs = coefs,
#                                   bval = t(x.normal[, !(colnames(x.normal) %in% reducing_coefs_inc)]))
#   pred_tumor_inc <- AgePred_tcga(name = paste0("reducing/tumor_inc_", perc, ".txt"),
#                                  coefs = coefs,
#                                  bval = t(x.tumor[, !(colnames(x.tumor) %in% reducing_coefs_inc)]))
#   print("after AgePred")
#   r_test_dec <- cor(y.test, pred_test_dec[1,])
#   r_normal_dec <- cor(y.normal, pred_normal_dec[1,])
#   r_tumor_dec <- cor(y.tumor, pred_tumor_dec[1,])
#   row_test_dec <- c(perc, r_test_dec, "Primary Tumour (Test)")
#   row_normal_dec <- c(perc, r_normal_dec, "Normal Tissue")
#   row_tumor_dec <- c(perc, r_tumor_dec, "Primary Tumour (Others)")
#   res_test_dec <- rbind(res_test_dec, row_test_dec)
#   res_normal_dec <- rbind(res_normal_dec, row_normal_dec)
#   res_tumor_dec <- rbind(res_tumor_dec, row_tumor_dec)
#   r_test_inc <- cor(y.test, pred_test_inc[1,])
#   r_normal_inc <- cor(y.normal, pred_normal_inc[1,])
#   r_tumor_inc <- cor(y.tumor, pred_tumor_inc[1,])
#   row_test_inc <- c(perc, r_test_inc, "Primary Tumour (Test)")
#   row_normal_inc <- c(perc, r_normal_inc, "Normal Tissue")
#   row_tumor_inc <- c(perc, r_tumor_inc, "Primary Tumour (Others)")
#   res_test_inc <- rbind(res_test_inc, row_test_inc)
#   res_normal_inc <- rbind(res_normal_inc, row_normal_inc)
#   res_tumor_inc <- rbind(res_tumor_inc, row_tumor_inc)
# }
# colnames(res_test_dec) <- c("Perc_missing", "R", "data")
# colnames(res_normal_dec) <- c("Perc_missing", "R", "data")
# colnames(res_tumor_dec) <- c("Perc_missing", "R", "data")
# write.csv(res_test_dec, paste0(dir_results, "model/test_reducing_dec.csv"))
# write.csv(res_normal_dec, paste0(dir_results, "model/normal_reducing_dec.csv"))
# write.csv(res_tumor_dec, paste0(dir_results, "model/tumor_reducing_dec.csv"))
# colnames(res_test_inc) <- c("Perc_missing", "R", "data")
# colnames(res_normal_inc) <- c("Perc_missing", "R", "data")
# colnames(res_tumor_inc) <- c("Perc_missing", "R", "data")
# write.csv(res_test_inc, paste0(dir_results, "model/test_reducing_inc.csv"))
# write.csv(res_normal_inc, paste0(dir_results, "model/normal_reducing_inc.csv"))
# write.csv(res_tumor_inc, paste0(dir_results, "model/tumor_reducing_inc.csv"))
# 
# # Plot results
# # # Load data
# # res_test <- read.csv("/Volumes/iorio/Irene/epiclock_dev/results/TCGA/model/test_reducing.csv")
# # res_normal <- read.csv("/Volumes/iorio/Irene/epiclock_dev/results/TCGA/model/normal_reducing.csv")
# # res_tumor <- read.csv("/Volumes/iorio/Irene/epiclock_dev/results/TCGA/model/tumor_reducing.csv")
# 
# # Plot
# color_type <- c("Primary Tumour (Test)" = "steelblue1",
#                 "Primary Tumour (Others)" = "steelblue3",
#                 "Normal Tissue" = "chartreuse3")
# res_dec <- as.data.frame(rbind(res_test_dec, res_normal_dec, res_tumor_dec))
# png(filename = "/group/iorio/Irene/epiclock/plots/reducing_coefficients_dec.png",
#     width = 12, height = 7, units = 'in', res = 600)
# ggplot(res_dec, aes(x = as.numeric(Perc_missing), y = as.numeric(R), fill = data, color = data)) +
#   geom_line() +
#   geom_point(shape=21, color="black", size=3) +
#   scale_fill_manual(values = color_type, name = "Data subset") +
#   scale_color_manual(values = color_type, name = "Data subset") +
#   theme_bw(base_size = 20) +
#   xlab("Missing Coefficients (%)") + ylab("Pearson Correlation (R)") 
# dev.off()
# 
# res_inc <- as.data.frame(rbind(res_test_inc, res_normal_inc, res_tumor_inc))
# png(filename = "/group/iorio/Irene/epiclock/plots/reducing_coefficients_inc.png",
#     width = 12, height = 7, units = 'in', res = 600)
# ggplot(res_inc, aes(x = as.numeric(Perc_missing), y = as.numeric(R), fill = data, color = data)) +
#   geom_line() +
#   geom_point(shape=21, color="black", size=3) +
#   scale_fill_manual(values = color_type, name = "Data subset") +
#   scale_color_manual(values = color_type, name = "Data subset") +
#   theme_bw(base_size = 20) +
#   xlab("Missing Coefficients (%)") + ylab("Pearson Correlation (R)") 
# dev.off()


# Performance by reducing the number of CpGs
res_test_dec <- c()
res_normal_dec <- c()
res_tumor_dec <- c()
res_test_inc <- c()
res_normal_inc <- c()
res_tumor_inc <- c()
reducing_perc <- seq(10, nrow(not0_coefs), by = 10) # reducing CpGs percentages
for(perc in reducing_perc){
  num_to_keep <- ceiling(perc*nrow(not0_coefs))
  # Coefficients to keep (the most extreme)
  reducing_coefs_dec <- not0_coefs_dec[seq_len(num_to_keep), 1]
  # Coefficients to keep (the least extreme)
  reducing_coefs_inc <- not0_coefs_inc[seq_len(num_to_keep), 1]
  pred_test_dec <- AgePred_tcga(name = paste0("reducing/test_dec_n", perc, ".txt"),
                                coefs = coefs,
                                bval = t(x.test[, !(colnames(x.test) %in% reducing_coefs_dec)]))
  pred_normal_dec <- AgePred_tcga(name = paste0("reducing/normal_dec_n", perc, ".txt"),
                                  coefs = coefs,
                                  bval = t(x.normal[, !(colnames(x.normal) %in% reducing_coefs_dec)]))
  pred_tumor_dec <- AgePred_tcga(name = paste0("reducing/tumor_dec_n", perc, ".txt"),
                                 coefs = coefs,
                                 bval = t(x.tumor[, !(colnames(x.tumor) %in% reducing_coefs_dec)]))
  pred_test_inc <- AgePred_tcga(name = paste0("reducing/test_inc_n", perc, ".txt"),
                                coefs = coefs,
                                bval = t(x.test[, !(colnames(x.test) %in% reducing_coefs_inc)]))
  pred_normal_inc <- AgePred_tcga(name = paste0("reducing/normal_inc_n", perc, ".txt"),
                                  coefs = coefs,
                                  bval = t(x.normal[, !(colnames(x.normal) %in% reducing_coefs_inc)]))
  pred_tumor_inc <- AgePred_tcga(name = paste0("reducing/tumor_inc_n", perc, ".txt"),
                                 coefs = coefs,
                                 bval = t(x.tumor[, !(colnames(x.tumor) %in% reducing_coefs_inc)]))
  print("after AgePred")
  r_test_dec <- cor(y.test, pred_test_dec[1,])
  r_normal_dec <- cor(y.normal, pred_normal_dec[1,])
  r_tumor_dec <- cor(y.tumor, pred_tumor_dec[1,])
  row_test_dec <- c(perc, r_test_dec, "Primary Tumour (Test)")
  row_normal_dec <- c(perc, r_normal_dec, "Normal Tissue")
  row_tumor_dec <- c(perc, r_tumor_dec, "Primary Tumour (Others)")
  res_test_dec <- rbind(res_test_dec, row_test_dec)
  res_normal_dec <- rbind(res_normal_dec, row_normal_dec)
  res_tumor_dec <- rbind(res_tumor_dec, row_tumor_dec)
  r_test_inc <- cor(y.test, pred_test_inc[1,])
  r_normal_inc <- cor(y.normal, pred_normal_inc[1,])
  r_tumor_inc <- cor(y.tumor, pred_tumor_inc[1,])
  row_test_inc <- c(perc, r_test_inc, "Primary Tumour (Test)")
  row_normal_inc <- c(perc, r_normal_inc, "Normal Tissue")
  row_tumor_inc <- c(perc, r_tumor_inc, "Primary Tumour (Others)")
  res_test_inc <- rbind(res_test_inc, row_test_inc)
  res_normal_inc <- rbind(res_normal_inc, row_normal_inc)
  res_tumor_inc <- rbind(res_tumor_inc, row_tumor_inc)
}
colnames(res_test_dec) <- c("Perc_missing", "R", "data")
colnames(res_normal_dec) <- c("Perc_missing", "R", "data")
colnames(res_tumor_dec) <- c("Perc_missing", "R", "data")
write.csv(res_test_dec, paste0(dir_results, "model/test_reducing_dec_n.csv"))
write.csv(res_normal_dec, paste0(dir_results, "model/normal_reducing_dec_n.csv"))
write.csv(res_tumor_dec, paste0(dir_results, "model/tumor_reducing_dec_n.csv"))
colnames(res_test_inc) <- c("Perc_missing", "R", "data")
colnames(res_normal_inc) <- c("Perc_missing", "R", "data")
colnames(res_tumor_inc) <- c("Perc_missing", "R", "data")
write.csv(res_test_inc, paste0(dir_results, "model/test_reducing_inc_n.csv"))
write.csv(res_normal_inc, paste0(dir_results, "model/normal_reducing_inc_n.csv"))
write.csv(res_tumor_inc, paste0(dir_results, "model/tumor_reducing_inc_n.csv"))

# Plot results
# # Load data
# res_test <- read.csv("/Volumes/iorio/Irene/epiclock_dev/results/TCGA/model/test_reducing.csv")
# res_normal <- read.csv("/Volumes/iorio/Irene/epiclock_dev/results/TCGA/model/normal_reducing.csv")
# res_tumor <- read.csv("/Volumes/iorio/Irene/epiclock_dev/results/TCGA/model/tumor_reducing.csv")

# Plot
color_type <- c("Primary Tumour (Test)" = "steelblue1",
                "Primary Tumour (Others)" = "steelblue3",
                "Normal Tissue" = "chartreuse3")
res_dec <- as.data.frame(rbind(res_test_dec, res_normal_dec, res_tumor_dec))
png(filename = "/group/iorio/Irene/epiclock/plots/reducing_coefficients_dec_n.png",
    width = 12, height = 7, units = 'in', res = 600)
ggplot(res_dec, aes(x = as.numeric(Perc_missing), y = as.numeric(R), fill = data, color = data)) +
  geom_line() +
  geom_point(shape=21, color="black", size=3) +
  scale_fill_manual(values = color_type, name = "Data subset") +
  scale_color_manual(values = color_type, name = "Data subset") +
  theme_bw(base_size = 20) +
  xlab("Number of Coefficients") + ylab("Pearson Correlation (R)") 
dev.off()

res_inc <- as.data.frame(rbind(res_test_inc, res_normal_inc, res_tumor_inc))
png(filename = "/group/iorio/Irene/epiclock/plots/reducing_coefficients_inc_n.png",
    width = 12, height = 7, units = 'in', res = 600)
ggplot(res_inc, aes(x = as.numeric(Perc_missing), y = as.numeric(R), fill = data, color = data)) +
  geom_line() +
  geom_point(shape=21, color="black", size=3) +
  scale_fill_manual(values = color_type, name = "Data subset") +
  scale_color_manual(values = color_type, name = "Data subset") +
  theme_bw(base_size = 20) +
  xlab("Number of Coefficients") + ylab("Pearson Correlation (R)") 
dev.off()

