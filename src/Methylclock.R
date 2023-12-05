# Methylclock

source("/group/iorio/Alessandro.D/EpiClock/src/dnaMethyAge.R")
library(methylclock)
library(dbplyr)

#chronological and biological DNAm age estimation
res<-DNAmAge(new_MethylationData, clocks=("Horvath"))
missCpGs <- checkClocks(new_MethylationData) # miss 25 CpGs from Horvarth model

# compute age acceleration
age1<-annotations_with_age2$Age
methylclock_age_acc<- DNAmAge(new_MethylationData, clocks=("Horvath"), age=age1, 
                         cell.count=F)
# plot correlation 
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




CMPcolnames <- as.data.frame(sort(colnames(filt_model_list)))






