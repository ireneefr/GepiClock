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























