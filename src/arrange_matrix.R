# arrange the right b-values df

colnames(MethylationData) [!colnames(MethylationData)== c("CPGs name","1331030", "906815", "1290907",
                                                          "908482", "905966", "905987") ]<- filt_model_list$model_name
a<-colnames(MethylationData)[startsWith(colnames(MethylationData), "1")]
b<-colnames(MethylationData)[startsWith(colnames(MethylationData), "9")]

# change also the names ofn the 6 uncommon samples
colnames(MethylationData)[colnames(MethylationData) %in% c("1331030", "906815", "1290907", "908482", "905966", "905987")] <- c("SC-1", "COLO-741", "HCC-56", "NCIH630", "SNB-19", "NCI-ADR-RES")
MethylationData2<-MethylationData

write.csv(filt_model_list, file= "../data/filt_model_list.csv", row.names=FALSE)
filt_model_list$model_name

# "1331030" "906815"  "1290907" "908482"  "905966"  "905987"  "1330932"
if (any(is.na(colnames(MethylationData)))) {
  cat("Ci sono valori mancanti nei nomi delle colonne.\n")
} else {
  cat("Non ci sono valori mancanti nei nomi delle colonne.\n")
}

cosm
tail(cosm)
head(cosm)
is.na(cosm)="TRUE"

if (any(is.na(cosm))) {
  cat("Ci sono valori mancanti nei nomi delle colonne.\n")
} else {
  cat("Non ci sono valori mancanti nei nomi delle colonne.\n")
}
num_na <- sum(is.na(cosm))



  
  
  
  
  