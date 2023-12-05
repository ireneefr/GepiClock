# Methylclock reference data test

library(tidyverse)
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










