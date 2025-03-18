####################################################
###   model age prediction CLs                   ###
###   Author: Digilio Alessandro                 ###
###   Date: 05/03/2025                           ###
####################################################

# in this script the performances of the new model vs the 4 most widely used clocks are evaluated and compared

library(tidyverse)
library(dplyr)

source("src/utils.R")
# source("src/load_data.R")
source("src/plot_generation_helpers.R")

set.seed(123)

# AGE PREDICTION
# compute age prediction with AgePred function (see: utils.R)
age_prediction <- AgePred_cell_lines(name = "cls_age_prediction.csv", coefs = coefficients_new_model, bval = b_values)

# check if the order of cls is the same
if(all(colnames(age_prediction) == cl_samples$CAccession)) {
  print("The order matches")
} else {
  print("There was an error in ordering. The orders do not match")
}

# add age prediction column to cl_samples information data frame
cl_samples <- cbind(cl_samples, age_prediction = age_prediction[1,])

# plot the distribution of the predicted age (from the model's CpGs) of the cell lines
plot_distribution_PredAge(cl_samples)

# compute median of the predicted age
median_age <- median(cl_samples$age_prediction, na.rm = TRUE)
# add age_group_median column
cl_samples$age_group_median <- ifelse(cl_samples$age_prediction > median_age, "Old", "Young")

