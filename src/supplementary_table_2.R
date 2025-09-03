#################################
###   Supplementary Table 2   ###
#################################

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Load model coefficients
model_coefs <- read.csv("results/TCGA/model/model_coefs.csv")
colnames(model_coefs) <- c("CpG", "Coefficient")

# Get annotation
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
model_coefs_ann <- merge(model_coefs[model_coefs$Coefficient != 0,], 
                         ann[,c("chr", "pos", "Name", "UCSC_RefGene_Name",
                               "UCSC_RefGene_Accession", "UCSC_RefGene_Group")],
                         by.x = "CpG", by.y = "Name", all.x = TRUE)
head(model_coefs_ann)

# Get Supplementary Table 2
writexl::write_xlsx(as.data.frame(model_coefs_ann), "Supplementary Tables/Supplementary_Table_2.xlsx")
