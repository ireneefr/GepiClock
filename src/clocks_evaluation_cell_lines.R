####################################################
###   Clocks Evaluation                          ###
###   Author: Digilio Alessandro                 ###
###   Date: 05/03/2025                           ###
####################################################

# in this script the performances of the new model vs the 4 most widely used clocks are evaluated and compared

library(dnaMethyAge)
library(tidyverse)
library(data.table)

src("src/plot_generation_helpers.R")

# load b-values and annotations from CMP of samples for which a"age at sampling" information is available
all_b_values_age_at_sampling <- fread("data/b_values_with_age_at_sampling/epiclock_CLs_methylationdata_.csv")
all_b_values_age_at_sampling <- all_b_values_age_at_sampling %>% column_to_rownames(var = "V1")
cl_samples_age_at_sampling <- read_csv("metadata/CMP_annotations/CMP_annotations_with_age_at_sampling/716_CLs_annotations.csv")
cl_samples_age_at_sampling <- cl_samples_age_at_sampling %>% column_to_rownames(var = "...1")
# check dimension 
dim(all_b_values_age_at_sampling) # 716 cell lines
dim(cl_samples_age_at_sampling) # 716 cell lines

# order in the same way
all_b_values_age_at_sampling <- all_b_values_age_at_sampling[, match(cl_samples_age_at_sampling$CAccession, colnames(all_b_values_age_at_sampling))]
# check if they have the same order
if (all(colnames(all_b_values_age_at_sampling) == cl_samples_age_at_sampling$CAccession)) {
  print("Cell lines are correctly aligned between b-values data and metadata.")
} else {
  stop("Error: Cell lines are not aligned.")
}

# compute the prediction with methyAge package 
availableClock() # list all supported clocks
# test the 4 widely used Epigenetic Clocks
Horvath <- "HorvathS2013"
Horvath2018 <- "HorvathS2018"
Levine <- "LevineM2018"
Hannum <- "HannumG2013"
horvath_age <- methyAge(as.matrix(all_b_values_age_at_sampling), clock=Horvath, inputation = FALSE, simple_mode = TRUE)
horvath2018_age <- methyAge(as.matrix(all_b_values_age_at_sampling), clock=Horvath2018, inputation = FALSE, simple_mode = TRUE)
levine_age <- methyAge(as.matrix(all_b_values_age_at_sampling), clock=Levine, inputation = FALSE, simple_mode = TRUE)
hannum_age <- methyAge(as.matrix(all_b_values_age_at_sampling), clock=Hannum, inputation = FALSE, simple_mode = TRUE)

# add information of age prediction of the 4 clocks
final_ann_ages <- cl_samples_age_at_sampling
final_ann_ages$Horvath2013 <- horvath_age$mAge
final_ann_ages$Horvath2018 <- horvath2018_age$mAge
final_ann_ages$Hannum <- hannum_age$mAge
final_ann_ages$Levine <- levine_age$mAge

# make them numeric
final_ann_ages$age_at_sampling <- as.numeric(final_ann_ages$age_at_sampling)
final_ann_ages$Horvath2013 <- as.numeric(final_ann_ages$Horvath2013)
final_ann_ages$Horvath2018 <- as.numeric(final_ann_ages$Horvath2018)
final_ann_ages$Hannum <- as.numeric(final_ann_ages$Hannum)
final_ann_ages$Levine <- as.numeric(final_ann_ages$Levine)

# PEARSON CORRELATION 
pearson_horvath <- cor(final_ann_ages$age_at_sampling, final_ann_ages$Horvath2013, method = "pearson", use = "complete.obs")
pearson_horvath2018 <- cor(final_ann_ages$age_at_sampling, final_ann_ages$Horvath2018, method = "pearson", use = "complete.obs")
pearson_hannum <- cor(final_ann_ages$age_at_sampling, final_ann_ages$Hannum, method = "pearson", use = "complete.obs")
pearson_levine <- cor(final_ann_ages$age_at_sampling, final_ann_ages$Levine, method = "pearson", use = "complete.obs")
# correlation df
clock_correlations <- data.frame(
  Clock = c("HorvathS2013", "HannumG2013", "HorvathS2018", "LevineM2018"),
  Pearson_Correlation = c(pearson_horvath, pearson_hannum, pearson_horvath2018, pearson_levine)
)
print(clock_correlations)

# plot clock performance correlations
plot_clocks_performances(clock_correlations)
