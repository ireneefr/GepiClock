# In this script there are various tests to try the dnaMethyAge and Methylclock predictors with
# the MethylationData obtained from the old preprocessing

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
                  "cg17274064", "cg21305265", "cg23180365", "cg24471894")# Methylclock

# METHYLCLOCK

# Methylclock reference data test

library(methylclock)
packageVersion("methylclock")
Ref_MethylationData<-read_csv("/group/iorio/Thanos/ReferenceData_MethyclockR/MethylationDataExample55.csv")
Ref_Prediction<-DNAmAge(Ref_MethylationData)
Ref_missCpGs <- checkClocks(Ref_MethylationData)

# 353 CpGs in Horvarth model
# Horvarth_CpGs<-read_csv("/group/iorio/Thanos/13059_2013_3156_MOESM3_ESM.csv")

# read annotation file (example)
covariates<-read_csv("/group/iorio/Thanos/ReferenceData_MethyclockR/SampleAnnotationExample55.csv")

# extract age vector and compute DNAmAge function to see the different to mAge and chronological age
age<-covariates$Age
age.example55 <- DNAmAge(Ref_MethylationData, age=age, 
                         cell.count=F)

meffil::meffil.list.cell.count.references()

?DNAmAge

age.example55
# ageAcc: Difference between DNAmAge and chronological age.
# ageAcc2: Residuals obtained after regressing chronological age and DNAmAge (similar to IEAA).
# ageAcc3: Residuals obtained after regressing chronological age and DNAmAge adjusted for cell counts (similar to EEAA)
# if you set cell.count = F you will have only ageAcc and ageAcc2

# Then, we can investigate, for instance, whether the accelerated age is associated with Autism. In that example we will
# use a non-parametric test (NOTE: use t-test or linear regression for large sample sizes)
autism <- covariates$diseaseStatus
?kruskal.test()
kruskal.test(age.example55$ageAcc.Horvath ~ autism)
kruskal.test(age.example55$ageAcc2.Horvath ~ autism)
# Kruskal-Wallis test: 
# Null Hypothesis (H0): The position parameters (means or medians) of the variable distributions are the same across all groups.
# Alternative Hypothesis (H1): At least one group has a position different from the others.
# interpretation of result: If the p-value associated with the test is sufficiently low (usually below a significance
# threshold, such as 0.05), one can reject the Null Hypothesis. This suggests that at least one group has a position 
# different from the others, indicating that there might be significant differences between the groups.
# (NOTE: use t-test or linear regression for large sample sizes)
# ex: look at glm model on reference tutorial...
# The association between disease status and DNAmAge estimated using Horvath’s method can be computed by:
# mod.skinHorvath <- glm (disease~ ageAcc2.Horvath , 
#                       data=age.example55,
#                       family="binomial")
# summary(mod.skinHorvath)


# Bayesian Neural Network to predict DNAmAge based on Horvath’s CpGs.
plotDNAmAge(age.example55$Horvath, age, tit="Bayesian Neural Network")
# plot correlation between model clock
plotCorClocks(age.example55)


dd <- GEOquery::getGEO("GSE19711")
gse19711 <- dd[[1]]
# The object gse19711is an ExpressionSet that can contains CpGs and phenotypic (e.g clinical) information
pheno <- pData(gse19711)
age <- as.numeric(pheno$`ageatrecruitment:ch1`)
disease <- pheno$`sample type:ch1`
disease
table(disease)

library(dbplyr)

#chronological and biological DNAm age estimation
res<-DNAmAge(new_MethylationData, clocks=("Horvath"))
missCpGs <- checkClocks(new_MethylationData) # miss 25 CpGs from Horvarth model

# compute age acceleration
age1<-annotations_with_age2$Age
methylclock_age_acc<- DNAmAge(new_MethylationData, clocks=("Horvath"), age=age1, 
                              cell.count=F)
# plot correlation 
plotDNAmAge(methylclock_age_acc$Horvath, age1)
# bayesian neural network
plotDNAmAge(methylclock_age_acc$Horvath, age1, tit="Bayesian Neural Network")

# MISSING CPGS:
# cg00374717 cg01873645 cg03167275 cg03565323 cg06121469 cg07337598 cg08251036 
# cg08434234 cg10523019 cg12768605 cg15661409 cg17099569 cg17655614 cg18440048
# cg22289837 cg25148589 cg27377450 cg01820374 cg04268405 cg11025793 cg14727952 
# cg17274064 cg21305265 cg23180365 cg24471894


# 
plot(res$id, horvath_age$Sample)
#
library(ggplot2)
ggplot(horvath_age, aes(x = mAge, y = res$Horvath)) +
  geom_point() +
  labs(title = "Methylclock vs dnaMethyAge predictions",
       x = "Methylclock mAge",
       y = "dnaMethyAge mAge")
# correlation test
cor.test(res$Horvath, horvath_age$mAge, method = c("pearson"))

# healthy samples from GSE68379 (with age annotation on CMP)
non_cancerous <- filt_model_list[filt_model_list$cancer_type == "Non-Cancerous",]
healthy_cosmic_ids <- non_cancerous$COSMIC_ID
healthy_cosmic_ids <- as.character(healthy_cosmic_ids)
q <- intersect(colnames(new_MethylationData), non_cancerous$COSMIC_ID)

#
healthy_MethylationData <- new_MethylationData[, q]
#
non_cancerous2<- non_cancerous[, c("COSMIC_ID", "tissue", "age_at_sampling", "cancer_type")]
names(non_cancerous2)[names(non_cancerous2) == "COSMIC_ID"] <- "Sample"
names(non_cancerous2)[names(non_cancerous2) == "age_at_sampling"] <- "Age"
# convert Age column as numeric
non_cancerous2[["Age"]] <- as.numeric(non_cancerous2[["Age"]])
non_cancerous2 <- non_cancerous2[-5,]

#
dim(healthy_MethylationData)
dim(non_cancerous2)
non_cancerous2 <- non_cancerous2[,-4]
non_cancerous_matrix<-as.matrix(non_cancerous2)



w <- as.data.frame(non_cancerous_matrix)
as.numeric(w$Age)
#
Horvarth <- "HorvathS2013"
horvath_age_healthy <- methyAge(healthy_MethylationData, clock=Horvarth, inputation = FALSE, simple_mode = TRUE)

# Predict epigenetic age and calculate age acceleration by using annotations from CMP (age_info parameter)
horvath_age_acc_healthy <- methyAge(healthy_MethylationData, clock=Horvarth, age_info=non_cancerous_matrix,
                                    inputation = FALSE, 
                                    simple_mode = TRUE,
                                    fit_method='Linear', do_plot=TRUE)
#
age_healthy<-non_cancerous2$Age
methylclock_age_acc_healthy<- DNAmAge(healthy_MethylationData, clocks=("Horvath"), age=age_healthy, 
                                      cell.count=F)
plotDNAmAge(methylclock_age_acc_healthy$Horvath, age_healthy, tit="Bayesian Neural Network")

# MODEL COMPARISON (blood samples)
# compare the Horvart, Levine and Hannum models for the blood samples with Methylclock package
blood_annotations <- subset(annotations_with_age2, tissue == "Haematopoietic and Lymphoid")

# extract blood_annotations cosmic_id and filter blood methylation data (b-values)
blood_cosmic_ids <- blood_annotations$Sample
blood_common_cosmic_ids <- intersect(colnames(new_MethylationData), blood_annotations$Sample)
blood_MethylationData <- new_MethylationData[, blood_common_cosmic_ids]
dim(blood_MethylationData)
dim(blood_annotations)

# run the prediction by making a comparison between the 3 models
blood_age1<-blood_annotations$Age
blood_methylclock_age_acc<- DNAmAge(blood_MethylationData, clocks=c("Horvath","Hannum", "Levine"), age=blood_age1, 
                                    cell.count=F)
plotCorClocks(blood_methylclock_age_acc)
