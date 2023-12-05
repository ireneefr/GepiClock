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