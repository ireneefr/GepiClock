# Build a df compatible with Methylclock package starting from beta values ("betas" matrix)

setwd("/group/iorio/Alessandro.D/EpiClock/pheb-master")



source("/group/iorio/Alessandro.D/EpiClock/src/data_preprocessing.R") # from here I have access to beta_values matrix (betas)
betas <-as.data.frame(betas)
betas <-cbind(rownames(betas),betas)
beta1<-t(betas)
beta1<-as.data.frame(beta1)

CPG1028celllines <-beta1[-c(1,2),]
str(CPG1028celllines)

# CPGs name have to be the first column
CPG1028celllines <- cbind(rownames(CPG1028celllines), CPG1028celllines)
rownames(CPG1028celllines)<-NULL
colnames(CPG1028celllines)[1] <- "CPGs name"


# save the MethylationData
library(tidyverse)
write.csv(CPG1028celllines, file= "../pheb-master/allCPG1028cellines.csv", row.names=FALSE)
MethylationData<-read_csv("../pheb-master/allCPG1028cellines.csv")
rownames(MethylationData)<-NULL

# remove the "X" in the colnames in order to have the same cosmicID of the annotation model list
col_names <- colnames(MethylationData)
new_col_names <- sub("^X", "", col_names)
colnames(MethylationData) <- new_col_names

# download updated annotations dmodel list
model_list<-read_csv("../data/model_list_20230923.csv")
# download cell lines information
celllines<-read.table("../data/celllines.txt", header = TRUE, sep = "\t")
# make a vector with all the cosmic id of celllines 
cosm<-as.character(celllines$COSMIC_ID) 

# filter the model_list only for the cell lines we want based on cosm
filt_model_list<-model_list[model_list$COSMIC_ID  %in% cosm, ]
dim(filt_model_list) # they should be 1028: this means that 6 are missing

# find the 6 missing cell lines
cosm <- cosm
annotationcosm<- filt_model_list$COSMIC_ID
common_ids <- cosm %in% annotationcosm
missing_ids <- cosm[!common_ids]
print("IDs missing in annotation cosm")
print(missing_ids)

a<- final1$COSMIC_ID
b<- csvPRED$id
common_ids <- b %in% a
missing_ids <- a[!common_ids]


# "1331030" "906815"  "1290907" "908482"  "905966"  "905987"  are missing
# "1290907" (HCC-56), "908482" (NCIH630), "905966" (SNB-19), "906815" (COLO-741), 
# "905987" (NCI-ADR-RES) are contaminated/misclassified

# Next step: check how many NAs in filt_model_list$age_at_sampling



# METHYLCLOCK PREDICTION
library(methylclock)
library(dbplyr)

# load the updated DNAmAge function from "scr"
source('/group/iorio/Alessandro.D/EpiClock/src/utils.R')

# test the clock with a small dataframe
# MethylationTest<-read_csv('/group/iorio/Alessandro.D/EpiClock/data/Methylation_test.csv')
# res<-DNAmAge(MethylationTest, clocks=("Horvath"))

#chronological and biological DNAm age estimation
res<-DNAmAge(MethylationData, clocks=("Horvath"))
res.f.cell.count<-DNAmAge(MethylationData, clocks=("Horvath"), cell.count=T
                          , cell.count.reference = "blood gse35069 complete")

# plot(res$Horvath,res.f.cell.count$Horvath)

missCpGs <- checkClocks(MethylationData)
write.csv(res, file= "../data/Horvarth_prediction.csv", row.names=FALSE)
res$id <- sub("^X", "", res$id)

######
#use<-model_list[model_list$COSMIC_ID %in% res$id,]

filt_new1<-filt_model_list[,c("age_at_sampling","model_name",
                              "COSMIC_ID","tissue")]
frid<-sum(is.na(filt_new1$age_at_sampling))

res1<-res
res1$id <- as.numeric(as.character(res1$id))
res1<-res1 %>% arrange(id) #from low to high id
res1<-res1[1:1022,]
str(res1) ####check if numeric 
dim(res1) #needs to be same dim as filt_new1
plot(res1$id)

filt_new1<-filt_new1 %>% arrange(COSMIC_ID)
plot(filt_new1$COSMIC_ID)

final<-cbind(filt_new1,res1)
str(final)


final1<-final[final$tissue=="Head and Neck",]
plot(final1$age_at_sampling,
     final1$Horvath,
     col=factor(final1$tissue),
     pch=16)

#install.packages("ggpubr")
library("ggpubr")
cor.test(final$age_at_sampling,
    final$Horvath, method = c("pearson", "kendall", "spearman"))
cor.test(x, y, method=c("pearson", "kendall", "spearman"))


renv::snapshot()

check<-function(data_frame) {
  any(is.na(colnames(data_frame))) }
result<-check(MethylationData)

# test number 2 
#restest<-res
#restest$id <- as.numeric(as.character(restest$id))
#restest<-restest %>% arrange(id)
#str(restest)
#plot(restest$id)
#id_test<-restest$id
#is_there_NA_test<-any(is.na(id_test))
#sum(is.na(id_test))
#write.table(restest, "../data/Horvart_Prediction.txt", row.names=F)

#try with numeric
id_MethylationDatatest<-colnames(MethylationData)
id_MethylationDatatest<-id_MethylationDatatest[-1]
id_MethylationDatatest <- as.numeric(as.character(id_MethylationDatatest))
str(id_MethylationDatatest)
str(id_test)
setdiff(id_MethylationDatatest,id_test)

#try with character
id_MethylationDatatest<-colnames(MethylationData)
id_MethylationDatatest<-id_MethylationDatatest[-1]
id_test<-res$id
str(id_MethylationDatatest)
str(id_test)
setdiff(id_MethylationDatatest,id_test)
?setdiff

#venn diagram to compare the 2 vectors
library(VennDiagram)
# Chart
a<-venn.diagram(
  x = list(id_test, id_MethylationDatatest),
  category.names = c("id from prediction" , "id from MethylationData"),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)
#
library(ggVennDiagram)
x <- list(
  A = id_MethylationDatatest, 
  B = id_test)
ggVennDiagram(x)


# Horvarth 353 cpgs
# colnames(Horvarth_CpGs) <- NULL
# colnames(Horvarth_CpGs) <- unname(Horvarth_CpGs[2,])
# Horvarth_CpGs<- Horvarth_CpGs[-c(1,2,3),]
Horvarth_CpGs_names <- Horvarth_CpGs$CpGmarker


















