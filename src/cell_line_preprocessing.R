# Cancer cell lines (GSE68379) pre-processing analysis
# Digilio Alessandro 

library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(ChAMP)
library(readxl)
library(tidyverse)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
RGset <- read.metharray.exp("/group/iorio/Alessandro.D/EpiClock/pheb-master/data/GSE68379_RAW/",
                            verbose=T) # read the test idat files from the test folder
# filter cross-reactive probes
xReactiveProbes <- read_excel("/group/iorio/Alessandro.D/EpiClock/pheb-master/metadata/48639-non-specific-probes-Illumina450k.xlsx") 

sentrix <- unlist(lapply(RGset@colData@rownames, function(x) paste(strsplit(x, "_")[[1]][-1],collapse="_")))
sample_sheet <- read_excel("/group/iorio/Alessandro.D/EpiClock/pheb-master/data/methSampleId_2_cosmicIds.xlsx")
sample_sheet$meth <- unlist(lapply(1:nrow(sample_sheet), function(x) paste(sample_sheet$Sentrix_ID[x],sample_sheet$Sentrix_Position[x], sep="_")))

cosmics <- unlist(lapply(sentrix, function(x) sample_sheet$cosmic_id[x == sample_sheet$meth]))
stopifnot(length(cosmics) == 1028)
all_cancer_type <- unique(sample_sheet$Tissue)
#
detP <- detectionP(RGset) # get detection pvalues
#
GRset.funnorm <- preprocessNoob(RGset, verbose = T) # "noob normalization"
Mset.funnorm <- mapToGenome(GRset.funnorm) # map to genome build
n_samples<-ncol(Mset.funnorm)
n_probes<-nrow(Mset.funnorm)
detP <- detP[match(featureNames(Mset.funnorm),rownames(detP)),] # ensure that probes are in same order than Mset

betas <- getBeta(Mset.funnorm)
NA_position <- detP > 0.05
betas[NA_position]<- NA

# remove cg with all NAs
# check thow many rowas we have with all NAs
rows_to_keep <- rowSums(is.na(betas)) / ncol(betas) < 1 # we don't have
betas <- betas[rows_to_keep,]
n_probes_keep <- nrow(betas)
# remove samples that fail >10% of probes
columns_to_keep <- colMeans(is.na(betas)) <= 0.1 
betas <- betas[, columns_to_keep]
n_samples_keep <- ncol(betas)
# remove cg with >20% NAs
rows_to_keep <- colMeans(is.na(betas)) <= 0.2 
betas <- betas[rows_to_keep,]
n_probes_keep2 <- nrow(betas)
#
keep <- !(rownames(betas) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")]) # filter out sex chromosomes
betas <- betas[keep,]
keep <- !(rownames(betas) %in% xReactiveProbes$TargetID)
betas <- betas[keep,]
n_probes_flt<-nrow(betas)

# imputation
n_nas<-sum(is.na(betas))
betas_imp<- champ.impute(betas, pd = NULL, method = "KNN")
rownames(myNorm)
myNorm <- champ.norm(betas_imp, arraytype = "450K", cores = 4, method = "BMIQ")
colnames(myNorm)<-sapply(strsplit(colnames(myNorm), split = "_"), function(x) x[1])

# save
rm(project_name)
rm(dir_data)
message("Saving betas preprocessed")
write.csv(myNorm, file= "/group/iorio/Alessandro.D/EpiClock/data/CLs_methylation_data.csv")
write(c(paste("Initial number of samples:", n_samples),
        paste("Samples with >10% failed probes:", n_samples-n_samples_keep),
        paste("Final number of samples:", n_samples_keep),
        paste("Initial number of probes:", n_probes),
        paste("Probes with all NAs:", n_probes-n_probes_keep),
        paste("Probes with >20% NAs:", n_probes_keep-n_probes_keep2),
        paste("Probes filtered:", n_probes_keep2-n_probes_flt),
        paste("Final number of probes:", n_probes_flt),
        paste("Number of NAs imputed:", n_nas)),
      file = paste("/group/iorio/Alessandro.D/EpiClock/data/info_preprocessing.txt", sep=""))


