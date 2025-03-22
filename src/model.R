#######################################################################
###   Elastic Net Regression to predict age in TCGA tumor samples   ###
###   Author: Irene Fernández Rebollo                               ###
###   Date: 25/02/2025                                              ###
#######################################################################

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
# Divide data in model (tumor), no model (tumor), and normal
samples_normal <- tcga_samples[tcga_samples$sample_type == "Solid Tissue Normal",]
samples_tumor <- tcga_samples[tcga_samples$sample_type %in% c("Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood"), ]
samples_model <- samples_tumor[!(samples_tumor$patient %in% c(samples_normal$patient, samples_tumor[duplicated(samples_tumor$patient), "patient"])),]
#samples_nomodel <- samples_tumor[!(samples_tumor$patient %in% samples_model$patient),]
bval_model <- tcga_bval[rownames(tcga_bval) %in% samples_model$barcode,]
# bval_nomodel <- tcga_bval[rownames(tcga_bval) %in% samples_nomodel$barcode,]
# bval_normal <- tcga_bval[rownames(tcga_bval) %in% samples_normal$barcode,]
toc()

tic("train-test")
# Train & Test split
set.seed(2025)
train_rows <- createDataPartition(samples_model[,"project_age"], p=0.8, list = FALSE)
train <- samples_model[train_rows, "barcode"]
test <- samples_model[-train_rows, "barcode"]

x.train <- bval_model[rownames(bval_model) %in% train, colnames(bval_model) %in% cpgs_shared]
x.test <- bval_model[rownames(bval_model) %in% test, colnames(bval_model) %in% cpgs_shared]
y.train <- samples_model[samples_model$barcode %in% train, "age_at_index"]
y.test <- samples_model[samples_model$barcode %in% test, "age_at_index"]
rm(tcga_bval, bval_model); gc() #free up memory

# Save barcodes in each set
writeLines(train, paste0(dir_results, "model/train.txt"))
writeLines(test, paste0(dir_results, "model/test.txt"))
# write.csv(x.train, paste0(dir_results, "model/x.train.csv"))
# write.csv(x.test, paste0(dir_results, "model/x.test.csv"))
# write.csv(y.train, paste0(dir_results, "model/y.train.csv"))
# write.csv(y.test, paste0(dir_results, "model/y.test.csv"))

toc()

tic("Create folds")
# Create stratified folds based on age and project
num_folds <- 10
folds <- createFolds(samples_model[train_rows, "project_age"], k = num_folds, list = FALSE)
toc()

tic("Elastic net")
# Alpha and lambda selection
fit <- cva.glmnet(x = x.train, y = y.train, foldid = folds, parallel = TRUE)
toc()

tic("Metrics")
# Obtain metrics 
print(fit$alpha)
df_alpha <- as.numeric(fit$alpha)
print(df_alpha)
df_lambda.1se <- c()
df_cvm <- c()
df_cvsd <- c()
df_nzero <- c()
for(i in 1:length(fit$alpha)){
  id <- which(fit$modlist[[i]]$lambda == fit$modlist[[i]]$lambda.1se)
  df_lambda.1se <- c(df_lambda.1se, fit$modlist[[i]]$lambda.1se)
  df_cvm <- c(df_cvm, fit$modlist[[i]]$cvm[id])
  df_cvsd <- c(df_cvsd, fit$modlist[[i]]$cvsd[id])
  df_nzero <- c(df_nzero, fit$modlist[[i]]$nzero[id])
}
df_fit <- data.frame(alpha = df_alpha,
                     lambda.1se = df_lambda.1se,
                     Measure = df_cvm,
                     SE = df_cvsd,
                     Nonzero = df_nzero)
write.csv(df_fit, paste0(dir_results, "model/stats_elasticnet.csv"))
toc()

tic("Elastic net with best alpha-lambda")
# Fit model with best alpha and lambda
best_id <- which(min(df_cvm) == df_cvm)
fit_best <- glmnet(x = x.train, y = y.train,
                   alpha = df_alpha[best_id],
                   lambda = df_lambda.1se[best_id])
saveRDS(fit_best, paste0(dir_results, "model/model.rds"))
model_coefs <- coef(fit_best)
write.csv(as.matrix(model_coefs), paste0(dir_results, "model/model_coefs.csv"))
toc()

stopCluster(cl)