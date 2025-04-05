##########################################################################
###   Compare Elastic Net trained with normal vs tumor (paired data)   ###
###   Author: Irene Fernández Rebollo                                  ###
###   Date: 25/02/2025                                                 ###
##########################################################################

setwd("/group/iorio/Irene/epiclock_dev")
dir_results <- "results/TCGA/"
source("src/utils.R")
library(doParallel)
library(parallel)
library(caret)
library(glmnet)
library(glmnetUtils)
library(tictoc)

# Parallelization
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Load data
tic("Load data")
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
toc()

tic("Divide data")
# Divide data in model_t (tumor paired), model_n (normal paired), and no model
samples_normal <- tcga_samples[tcga_samples$sample_type == "Solid Tissue Normal",]
samples_tumor <- tcga_samples[tcga_samples$sample_type == "Primary Tumor", ]
samples_tumor <- samples_tumor[!duplicated(samples_tumor$patient),] #remove duplicates
patients_paired <- intersect(samples_normal$patient, samples_tumor$patient)
samples_model_n <- samples_normal[samples_normal$patient %in% patients_paired,]
samples_model_t <- samples_tumor[samples_tumor$patient %in% patients_paired,]
#samples_nomodel <- rbind(samples_tumor[!(samples_tumor$patient %in% samples_model_t$patient),],
#                         tcga_samples[tcga_samples$sample_type == "Primary Blood Derived Cancer - Peripheral Blood", ])
bval_model_t <- tcga_bval[rownames(tcga_bval) %in% samples_model_t$barcode,]
bval_model_n <- tcga_bval[rownames(tcga_bval) %in% samples_model_n$barcode,]
#bval_nomodel <- tcga_bval[rownames(tcga_bval) %in% samples_nomodel$barcode,]
toc()

# Train & Test split
tic("train-test")
set.seed(2025)

train_rows_t <- createDataPartition(samples_model_t[,"project_age"], p=0.8, list = FALSE)
train_t <- samples_model_t[train_rows_t, "barcode"]
test_t <- samples_model_t[-train_rows_t, "barcode"]

train_rows_n <- createDataPartition(samples_model_n[,"project_age"], p=0.8, list = FALSE)
train_n <- samples_model_n[train_rows_n, "barcode"]
test_n <- samples_model_n[-train_rows_n, "barcode"]

x.train_t <- bval_model_t[rownames(bval_model_t) %in% train_t, colnames(bval_model_t) %in% cpgs_shared]
x.test_t <- bval_model_t[rownames(bval_model_t) %in% test_t, colnames(bval_model_t) %in% cpgs_shared]
y.train_t <- samples_model_t[samples_model_t$barcode %in% train_t, "age_at_index"]
y.test_t <- samples_model_t[samples_model_t$barcode %in% test_t, "age_at_index"]

x.train_n <- bval_model_n[rownames(bval_model_n) %in% train_n, colnames(bval_model_n) %in% cpgs_shared]
x.test_n <- bval_model_n[rownames(bval_model_n) %in% test_n, colnames(bval_model_n) %in% cpgs_shared]
y.train_n <- samples_model_n[samples_model_n$barcode %in% train_n, "age_at_index"]
y.test_n <- samples_model_n[samples_model_n$barcode %in% test_n, "age_at_index"]
rm(tcga_bval, bval_model_t, bval_model_n); gc() #free up memory

# Save barcodes in each set
writeLines(train_t, paste0(dir_results, "model/train_t.txt"))
writeLines(test_t, paste0(dir_results, "model/test_t.txt"))
writeLines(train_n, paste0(dir_results, "model/train_n.txt"))
writeLines(test_n, paste0(dir_results, "model/test_n.txt"))
# write.csv(x.train_t, paste0(dir_results, "model/x.train_t.csv"))
# write.csv(x.test_t, paste0(dir_results, "model/x.test_t.csv"))
# write.csv(y.train_t, paste0(dir_results, "model/y.train_t.csv"))
# write.csv(y.test_t, paste0(dir_results, "model/y.test_t.csv"))
# 
# write.csv(x.train_n, paste0(dir_results, "model/x.train_n.csv"))
# write.csv(x.test_n, paste0(dir_results, "model/x.test_n.csv"))
# write.csv(y.train_n, paste0(dir_results, "model/y.train_n.csv"))
# write.csv(y.test_n, paste0(dir_results, "model/y.test_n.csv"))
toc()

tic("Create folds")
# Create stratified folds based on age and project
num_folds <- 10
folds_t <- createFolds(samples_model_t[train_rows_t, "project_age"], k = num_folds, list = FALSE)
folds_n <- createFolds(samples_model_n[train_rows_n, "project_age"], k = num_folds, list = FALSE)
toc()

tic("Elastic net")
# Alpha and lambda selection
fit_t <- cva.glmnet(x = x.train_t, y = y.train_t, foldid = folds_t, parallel = TRUE)
fit_n <- cva.glmnet(x = x.train_n, y = y.train_n, foldid = folds_n, parallel = TRUE)
toc()


tic("Metrics")
# Obtain metrics 
df_alpha_t <- as.numeric(fit_t$alpha)
df_lambda.1se_t <- c()
df_cvm_t <- c()
df_cvsd_t <- c()
df_nzero_t <- c()
for(i in 1:length(fit_t$alpha)){
  id_t <- which(fit_t$modlist[[i]]$lambda == fit_t$modlist[[i]]$lambda.1se)
  df_lambda.1se_t <- c(df_lambda.1se_t, fit_t$modlist[[i]]$lambda.1se)
  df_cvm_t <- c(df_cvm_t, fit_t$modlist[[i]]$cvm[id_t])
  df_cvsd_t <- c(df_cvsd_t, fit_t$modlist[[i]]$cvsd[id_t])
  df_nzero_t <- c(df_nzero_t, fit_t$modlist[[i]]$nzero[id_t])
}
df_fit_t <- data.frame(alpha = df_alpha_t,
                       lambda.1se = df_lambda.1se_t,
                       Measure = df_cvm_t,
                       SE = df_cvsd_t,
                       Nonzero = df_nzero_t
)
write.csv(df_fit_t, paste0(dir_results, "model/stats_elasticnet_t.csv"))

df_alpha_n <- as.numeric(fit_n$alpha)
df_lambda.1se_n <- c()
df_cvm_n <- c()
df_cvsd_n <- c()
df_nzero_n <- c()
for(i in 1:length(fit_n$alpha)){
  id_n <- which(fit_n$modlist[[i]]$lambda == fit_n$modlist[[i]]$lambda.1se)
  df_lambda.1se_n <- c(df_lambda.1se_n, fit_n$modlist[[i]]$lambda.1se)
  df_cvm_n <- c(df_cvm_n, fit_n$modlist[[i]]$cvm[id_n])
  df_cvsd_n <- c(df_cvsd_n, fit_n$modlist[[i]]$cvsd[id_n])
  df_nzero_n <- c(df_nzero_n, fit_n$modlist[[i]]$nzero[id_n])
}
df_fit_n <- data.frame(alpha = df_alpha_n,
                       lambda.1se = df_lambda.1se_n,
                       Measure = df_cvm_n,
                       SE = df_cvsd_n,
                       Nonzero = df_nzero_n
)
write.csv(df_fit_n, paste0(dir_results, "model/stats_elasticnet_n.csv"))
toc()

tic("Elastic net with best alpha-lambda")
# Fit model with best alpha and lambda
best_id_t <- which(min(df_cvm_t) == df_cvm_t)
best_id_n <- which(min(df_cvm_n) == df_cvm_n)
fit_best_t <- glmnet(x = x.train_t, y = y.train_t,
                     alpha = df_alpha_t[best_id_t],
                     lambda = df_lambda.1se_t[best_id_t])
fit_best_n <- glmnet(x = x.train_n, y = y.train_n,
                     alpha = df_alpha_n[best_id_n],
                     lambda = df_lambda.1se_n[best_id_n])
saveRDS(fit_best_t, file = paste0(dir_results, "model/model_t.rds"))
saveRDS(fit_best_n, file = paste0(dir_results, "model/model_n.rds"))
model_coefs_t <- coef(fit_best_t)
model_coefs_n <- coef(fit_best_n)
write.csv(as.matrix(model_coefs_t), paste0(dir_results, "model/model_coefs_t.csv"))
write.csv(as.matrix(model_coefs_t), paste0(dir_results, "model/model_coefs_n.csv"))
toc()

stopCluster(cl)

