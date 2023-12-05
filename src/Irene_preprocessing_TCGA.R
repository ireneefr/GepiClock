############################################
###   Preprocess TCGA methylation data   ###
###   Author: Irene Fernández Rebollo    ###
###   Date: 30/10/2023                   ###
############################################

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(SummarizedExperiment)
library(ChAMP)

preprocess_tcga <- function(project_name, dir_data){
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  xReactiveProbes <- read.csv("metadata/48639-non-specific-probes-Illumina450k.csv")
  print(paste0("Loading: ", paste(dir_data, project_name, "/data_met_raw.RData", sep="")))
  load(paste(dir_data, project_name, "/data_met_raw.RData", sep=""))
  betas <- as.data.frame(assay(data_met))
  # remove cg with all NAs
  rows_to_keep <- rowSums(is.na(betas)) / ncol(betas) < 1 
  betas <- betas[rows_to_keep,]
  # remove samples that fail >10% of probes
  columns_to_keep <- colMeans(is.na(betas)) <= 0.1 
  betas <- betas[, columns_to_keep]
  # remove probes: crossreactive, in sexual chromosomes and not in ann450k (rs)
  probes_flt <- c(xReactiveProbes$TargetID,
                  ann450k$Name[ann450k$chr %in% c("chrX","chrY")],
                  rownames(betas)[!(rownames(betas) %in% ann450k$Name)])
  betas_flt <- betas[!(rownames(betas) %in% probes_flt),]
  # only keep cg probes
  betas_flt_cg <- betas_flt[grep("^cg", rownames(betas_flt)),]
  # imputation
  betas_flt_cg$mean <- rowMeans(as.matrix(betas_flt_cg), na.rm = TRUE)
  nas <- which(is.na(betas_flt_cg), arr.ind = TRUE)
  for(i in 1:nrow(nas)) {
    betas_flt_cg[nas[i, "row"], nas[i, "col"]] <- betas_flt_cg[nas[i, "row"], "mean"]
  }
  betas_flt_cg$mean <- NULL
  # BMIQ
  betas_bmiq <- champ.norm(betas_flt_cg, arraytype = "450K", cores = 1, method = "BMIQ")
  # Save
  save(betas_bmiq, file = paste(dir_data, project_name, "/data_met_preprocessed.RData", sep=""))
}