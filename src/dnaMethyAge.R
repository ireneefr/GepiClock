# dnaMethyAge 

library(tidyverse)
library(dnaMethyAge)

# import Methylation data
MethylationData<-read_csv("../pheb-master/allCPG1028cellines.csv")
rownames(MethylationData)<-NULL

# remove the "X" in the colnames in order to have the same cosmicID of CMP
col_names <- colnames(MethylationData)
new_col_names <- sub("^X", "", col_names)
colnames(MethylationData) <- new_col_names

# CpGs should be the rownames of the df and MethylationData has to be a matrix
rownames(MethylationData) <- MethylationData$`CPGs name`
rownames_original <- rownames(MethylationData)
MethylationData2 <- MethylationData[, -1]
rownames(MethylationData2) <- rownames_original
MethylationData2 <- as.matrix(MethylationData2)

# import annotations from CMP (filter only the ones with annotated age at sampling)
annotations_with_age <- read.table("/group/iorio/Alessandro.D/EpiClock/data/annotations_with_age_annotated.txt", header = TRUE)
annotations_with_age2 <- annotations_with_age[, c("COSMIC_ID", "tissue", "age_at_sampling")]

# change the colnames in order to fit with MethyAge function
names(annotations_with_age2)[names(annotations_with_age2) == "COSMIC_ID"] <- "Sample"
names(annotations_with_age2)[names(annotations_with_age2) == "age_at_sampling"] <- "Age"
# convert Age column as numeric
annotations_with_age2[["Age"]] <- as.numeric(annotations_with_age2[["Age"]])

# filter b-values of samples that has annotated age in CMP
common_cosmic_ids <- intersect(colnames(MethylationData2), annotations_with_age$COSMIC_ID)
new_MethylationData <- MethylationData2[, common_cosmic_ids]
new_MethylationData <- as.matrix(new_MethylationData)

# are there NAs?
any(is.na(annotations_with_age2$Age)) #TRUE 
na_position <- which(is.na(annotations_with_age2), arr.ind = TRUE)
print(na_position)
cosmic_NAs <- annotations_with_age2$Sample[na_position]
print(cosmic_NAs) 
#  910944 1330933  753620  907170  908440  908447  910567  910941  949172  949174 1240181  907298  683667  907299  907050
# 907050  907050  907050  907050  907050  907050  907050  907050  907050  907050  907050  907050  907050 are the NAs sample
# remove NAs from both the matrixes
NAs_samples <- c(
  "910944", "1330933", "753620", "907170", "908440",
  "908447", "910567", "910941", "949172", "949174",
  "1240181", "907298", "683667", "907299", "907050")
# remove from b-values matrix
new_MethylationData <- new_MethylationData[, !colnames(new_MethylationData) %in% NAs_samples]
# remove from annotation matrix
NAs_index <- which(annotations_with_age2$Sample %in% NAs_samples)
annotations_with_age2 <- annotations_with_age2[-NAs_index,]
annotations_with_age2<- annotations_with_age2[!(annotations_with_age2$Sample %in% NAs_index),]


# check dimension
# check dimension in order to have the same samples
dim(new_MethylationData) # 730 cell lines
dim(annotations_with_age2) # 730 cell lines

# compute the prediction
availableClock() # list all supported clocks
Horvath <- "HorvathS2013"
horvath_age <- methyAge(new_MethylationData, clock=Horvath, inputation = FALSE, simple_mode = TRUE)

# Predict epigenetic age and calculate age acceleration by using annotations from CMP (age_info parameter)
horvath_age_acc <- methyAge(new_MethylationData, clock=Horvath, age_info=annotations_with_age2, inputation = FALSE, 
                            simple_mode = TRUE,
                            fit_method='Linear', do_plot=TRUE)

# MISSING CPGS:
# cg00374717 cg01873645 cg03167275 cg03565323 cg06121469 cg07337598 cg08251036 
# cg08434234 cg10523019 cg12768605 cg15661409 cg17099569 cg17655614 cg18440048
# cg22289837 cg25148589 cg27377450 cg01820374 cg04268405 cg11025793 cg14727952 
# cg17274064 cg21305265 cg23180365 cg24471894
missing_CpGs <- c("cg00374717", "cg01873645", "cg03167275", "cg03565323", "cg06121469", "cg07337598", "cg08251036", 
                  "cg08434234", "cg10523019", "cg12768605", "cg15661409", "cg17099569", "cg17655614", "cg18440048",
                  "cg22289837", "cg25148589", "cg27377450", "cg01820374", "cg04268405", "cg11025793", "cg14727952", 
                  "cg17274064", "cg21305265", "cg23180365", "cg24471894")





