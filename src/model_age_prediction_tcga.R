#####################################################
###   Prediction epigenetic clocks in TCGA data   ###
###   Author: Irene Fernández Rebollo             ###
###   Date: 25/02/2025                            ###
#####################################################

source("/group/iorio/Irene/git_epiclock/src/utils.R")
BiocManager::install("dnaMethyAge")
library(dnaMethyAge)
dir_data = "/group/iorio/Irene/epiclock/data/"

# Laod data
projects_tcga <- TCGA_projects()
tcga_bval <- TCGA_Bvalues()

# HorvathS2013 predictions
horvath_pred <- methyAge(t(tcga_bval), clock="HorvathS2013", inputation = FALSE, simple_mode = TRUE)
write.csv(horvath_pred, "/group/iorio/Irene/git_epiclock/res/predictions/pred_HorvathS2013.csv")

# HannumG2013 predictions
hannum_pred <- methyAge(t(tcga_bval), clock="HannumG2013", inputation = FALSE, simple_mode = TRUE)
write.csv(hannum_pred, "/group/iorio/Irene/git_epiclock/res/predictions/pred_HannumG2013.csv")

# HorvathS2018 predictions
horvath2_pred <- methyAge(t(tcga_bval), clock="HorvathS2018", inputation = FALSE, simple_mode = TRUE)
write.csv(horvath2_pred, "/group/iorio/Irene/git_epiclock/res/predictions/pred_HorvathS2018.csv")

# LevineM2018 predictions
levine_pred <- methyAge(t(tcga_bval), clock="LevineM2018", inputation = FALSE, simple_mode = TRUE)
write.csv(levine_pred, "/group/iorio/Irene/git_epiclock/res/predictions/pred_LevineM2018.csv")

# RebolloI2025 predictions
coefs <- read.csv("/group/iorio/Irene/git_epiclock/res/model/model_coefs.csv")
AgePred_tcga(name = "pred_RebolloI2025.csv",
        coefs = coefs,
        bval = t(tcga_bval))

# Model paired tumor predictions
coefs <- read.csv("/group/iorio/Irene/git_epiclock/res/model/model_coefs_t.csv")
AgePred_tcga(name = "pred_paired_t.csv",
        coefs = coefs,
        bval = t(tcga_bval))

# Model paired normal predictions
coefs <- read.csv("/group/iorio/Irene/git_epiclock/res/model/model_coefs_n.csv")
AgePred_tcga(name = "pred_paired_n.csv",
        coefs = coefs,
        bval = t(tcga_bval))
