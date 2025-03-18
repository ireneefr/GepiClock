####################################################
###   Utils                                      ###
###   Author: Digilio Alessandro                 ###
###   Date: 05/03/2025                           ###
####################################################

# in this script are collected all the functions used
library(data.table)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(fgsea)
library(writexl)

source("src/global_params.R")

############ load_model_coefficients() ###################
# the 'load_model_coefficients' function automatically load and filter model coefficients
load_model_coefficients <- function(file_path) {
  coefficients <- read_csv(file_path)
  
  new_model <- coefficients[coefficients$s0 != 0, ]
  new_model <- new_model[-1, ]
  model_CpGs <- new_model$...1
  model_CpGs <- as.character(model_CpGs)
  
  coefficients_new_model <- as.data.frame(new_model)
  row.names(coefficients_new_model) <- coefficients_new_model$...1
  
  coefficients_new_model$...1 <- NULL
  
  return(list(
    model_CpGs = model_CpGs,
    coefficients_new_model = coefficients_new_model
  ))
}





############ load_b_values_and_annotations() ###################
# the 'load_b_values_and_annotations' function automatically loads 
# cell lines' methylation data (b-values) and annotations, and adds 
# doubling time from experimental settings.

load_b_values_and_annotations <- function(b_values_path, model_CpGs, iorio_path, annotations_path, experimental_settings_path) { 
  
  # ============================
  # LOAD METHYLATION DATA (B-VALUES)
  # ============================
  
  all_b_values <- fread(b_values_path) %>% column_to_rownames(var = "V1")
  b_values <- all_b_values[rownames(all_b_values) %in% model_CpGs, , drop = FALSE] # filter only model CpGs
  rownames(b_values) <- model_CpGs
  
  # ============================
  # LOAD CELL LINE ANNOTATIONS
  # ============================
  
  # load Iorio et al. (2016) cell line annotations
  iorio <- read.table(file = iorio_path, header = TRUE, sep = "\t")
  # load all model list CMP annotations
  cl_samples <- read.csv(annotations_path)
  cl_samples <- cl_samples %>% filter(model_name %in% iorio$Title.2) # filter only cell lines in Iorio dataset
  cl_samples$CAccession <- iorio$CAccession[match(cl_samples$model_name, iorio$Title.2)]
  # add acronyms for cancer types
  cl_samples$acronyms_cancer_type <- unlist(lapply(cl_samples$cancer_type, function(x) acronym_map[[x]]))
  # change model_id column name to match GDSC dataset
  colnames(cl_samples)[colnames(cl_samples) == "model_id"] <- "SANGER_MODEL_ID"
  
  # ============================
  # LOAD EXPERIMENTAL SETTINGS (DOUBLING TIME)
  # ============================
  
  # load experimental growth rate data
  all_cl_experimental_settings <- read.table(experimental_settings_path, header = TRUE)
  # rename column for matching
  colnames(all_cl_experimental_settings)[colnames(all_cl_experimental_settings) == "model_id"] <- "SANGER_MODEL_ID"
  # filter only relevant cell lines
  cl_experimental_settings <- all_cl_experimental_settings %>% 
    filter(SANGER_MODEL_ID %in% cl_samples$SANGER_MODEL_ID) %>%
    select(SANGER_MODEL_ID, doubling_time_hours) # keep only doubling time column
  # merge doubling time into cl_samples
  cl_samples <- left_join(cl_samples, cl_experimental_settings, by = "SANGER_MODEL_ID")
  
  # ============================
  # FILTER B-VALUES TO MATCH CELL LINES
  # ============================
  
  b_values <- b_values[, colnames(b_values) %in% cl_samples$CAccession, drop = FALSE] # keep only matching cell lines
  cl_samples <- cl_samples %>% filter(!duplicated(CAccession)) 
  b_values <- b_values[, match(cl_samples$CAccession, colnames(b_values))]
  
  # ensure the order matches
  if (all(colnames(b_values) == cl_samples$CAccession)) {
    message("The order of cell lines in b_values and cl_samples matches.")
  } else {
    stop("Error: The orders of b_values and cl_samples do not match.")
  }
  
  # ============================
  # RETURN PROCESSED DATA
  # ============================
  
  return(list(
    b_values = b_values,
    cl_samples = cl_samples
  ))
}





############ load_450k_methylation_data() ###################
# the the 'load_450k_methylation_data' function load and process 450k CpG methylation data (b-values) with cell line annotations
# this function loads the cell line methylation dataset (b-values), filters it to retain only cell lines present in 
# the Iorio et al. (2016) dataset, and merges experimental doubling time data. The function ensures that 
# the methylation data and cell line metadata are correctly aligned.this function it's used only in UMAP_methylation_data to compare 
# the methylation profile pattenrs with the ones related to the new's model CpGS
load_450k_methylation_data <- function(all_b_values_path, iorio_path, annotations_path, experimental_settings_path) {
  
  # ============================
  # LOAD METHYLATION DATA (B-VALUES)
  # ============================
  all_b_values <- fread(all_b_values_path) %>% column_to_rownames(var = "V1")
  
  # ============================
  # LOAD CELL LINE ANNOTATIONS
  # ============================
  
  # load Iorio et al. (2016) cell line annotations
  iorio <- read.table(file = iorio_path, header = TRUE, sep = "\t")
  
  # load all model list CMP annotations
  cl_samples <- read.csv(annotations_path)
  
  # filter only cell lines present in the Iorio dataset
  cl_samples <- cl_samples %>% filter(model_name %in% iorio$Title.2)
  
  # match and assign CAccession from Iorio dataset
  cl_samples$CAccession <- iorio$CAccession[match(cl_samples$model_name, iorio$Title.2)]
  
  # add acronyms for cancer types (using a predefined mapping list `acronym_map`)
  cl_samples$acronyms_cancer_type <- unlist(lapply(cl_samples$cancer_type, function(x) acronym_map[[x]]))
  
  # rename model_id column to match GDSC dataset
  colnames(cl_samples)[colnames(cl_samples) == "model_id"] <- "SANGER_MODEL_ID"
  
  # ============================
  # LOAD EXPERIMENTAL SETTINGS (DOUBLING TIME)
  # ============================
  
  # load experimental growth rate data
  all_cl_experimental_settings <- read.table(experimental_settings_path, header = TRUE)
  
  # rename model_id column for consistency
  colnames(all_cl_experimental_settings)[colnames(all_cl_experimental_settings) == "model_id"] <- "SANGER_MODEL_ID"
  
  # filter only relevant cell lines
  cl_experimental_settings <- all_cl_experimental_settings %>%
    filter(SANGER_MODEL_ID %in% cl_samples$SANGER_MODEL_ID) %>%
    select(SANGER_MODEL_ID, doubling_time_hours) # Keep only doubling time column
  
  # merge doubling time into cl_samples metadata
  cl_samples <- left_join(cl_samples, cl_experimental_settings, by = "SANGER_MODEL_ID")
  
  # ============================
  # FILTER B-VALUES TO MATCH CELL LINES
  # ============================
  
  # keep only matching cell lines in methylation data
  all_b_values <- all_b_values[, colnames(all_b_values) %in% cl_samples$CAccession, drop = FALSE]
  
  # remove duplicates in metadata
  cl_samples <- cl_samples %>% filter(!duplicated(CAccession))
  
  # ensure matching order between methylation matrix and metadata
  all_b_values <- all_b_values[, match(cl_samples$CAccession, colnames(all_b_values))]
  
  # verify that the order is correctly aligned
  if (all(colnames(all_b_values) == cl_samples$CAccession)) {
    message("he order of cell lines in all_b_values and cl_samples matches.")
  } else {
    stop("Error: The orders of all_b_values and cl_samples do not match.")
  }
  
  # ============================
  # RETURN PROCESSED DATA
  # ============================
  
  return(list(
    all_b_values = all_b_values,
    cl_samples = cl_samples
  ))
}





############ AgePred() ###################
# the 'AgePred' function predicts biological age of cell lines using CpG methylation values based on a model’s coefficients
AgePred_cell_lines <- function(dir_predictions = "/group/iorio/Alessandro.D/EpiClock/results/predictions/",
                    name, # file name
                    coefs, # coefficient of the model's CpGs
                    bval){ # b_values
  
  not0_coefs <- rownames(coefs)[coefs[,1] != 0]
  cpgs_common <- intersect(not0_coefs, rownames(bval))
  cpgs_missing <- not0_coefs[!(not0_coefs %in% cpgs_common)][-1]
  message(paste("There are", length(cpgs_missing), "CpGs missing:",
                paste(cpgs_missing, sep = ", ")))
  
  pred <- t(coefs[cpgs_common, 1]) %*% as.matrix(bval[cpgs_common,])
  pred <- coefs[1, 1] + pred
  write.csv(pred, paste0(dir_predictions, name))
  
  return(pred)
}





############ filter_tissue_or_cancer() ###################
# the'filter_tissue_or_cancer' function filters a dataset based on cancer type or tissue availability, selecting
# categories with a specified minimum number of cell lines and classifying each entry by its cancer type or tissue
filter_tissue_or_cancer <- function(df, cancer_col, tissue_col, accession_col, threshold) {
  new_col <- rep(NA, nrow(df)) # empty column
  cancer_types_to_keep <- tapply(df[[accession_col]], INDEX = as.factor(df[[cancer_col]]), 
                                 function(i) length(unique(i))) >= threshold # check each cancer type to see if there are more than 'threshold' cell lines
  
  pos_to_confirm <- which(df[[cancer_col]] %in% names(cancer_types_to_keep)[cancer_types_to_keep])
  new_col[pos_to_confirm] <- df[[cancer_col]][pos_to_confirm] # populate the new column with cancer types that meet the threshold
  few_cancer_df <- df[is.na(new_col), ]  # filter records with cancer types under the threshold
  tissue_types_to_keep <- tapply(few_cancer_df[[accession_col]], INDEX = as.factor(few_cancer_df[[tissue_col]]), 
                                 function(i) length(unique(i))) >= threshold  # check each tissue type to see if there are more than 'threshold' cell lines
  
  pos_to_confirm_tissue <- which(df[[tissue_col]] %in% names(tissue_types_to_keep)[tissue_types_to_keep] & is.na(new_col))
  new_col[pos_to_confirm_tissue] <- df[[tissue_col]][pos_to_confirm_tissue]  # populate the new column with the tissue name if it meets the threshold and isn't already populated
  
  df$tissue_or_cancer_chosen <- new_col   # add the new column to the input dataframe
  
  return(df)
}






############ generate_merged_drug_dataset() ###################
# the'generate_merged_drug_dataset' function filters a dataset based on cancer type or tissue availability, selecting
# categories with a specified minimum number of cell lines and classifying each entry by its cancer type or tissue
generate_merged_drug_dataset <- function(GDSC_combined, cl_samples, threshold = 10) {
  # merge GDSC_combined with relevant columns from cl_samples
  GDSC_age <- left_join(
    GDSC_combined, 
    cl_samples[, c("SANGER_MODEL_ID", "age_prediction", "age_group_median", "age_at_sampling",
                   "tissue", "msi_status", "cancer_type","acronyms_cancer_type", "CAccession")],
    by = "SANGER_MODEL_ID"
  ) %>%
    filter(!is.na(age_prediction) & !is.na(tissue) & !is.na(msi_status))
  
  # apply filter_tissue_or_cancer function to classify based on cancer type or tissue
  GDSC_age <- filter_tissue_or_cancer(GDSC_age, "cancer_type", "tissue", "CAccession", threshold)
  
  # update missing values in tissue_or_cancer_chosen column
  GDSC_age$tissue_or_cancer_chosen[is.na(GDSC_age$tissue_or_cancer_chosen)] <- "others"
  
  # create status column to classify entries as Cancer Type, Tissue Type, or Others
  GDSC_age$status <- ifelse(GDSC_age$tissue_or_cancer_chosen == "others", "Others",
                            ifelse(GDSC_age$tissue_or_cancer_chosen %in% GDSC_age$cancer_type, "Cancer Type", "Tissue Type"))
  
  # create column for cancer specific analysys
  GDSC_age$cancer_specific <- "Excluded"
  
  # mark as "Included as Cancer Type" where status is "Cancer Type" in order to include onl ythe sample with > 10 cancer types for the cancer type specific analysis
  GDSC_age$cancer_specific[GDSC_age$status == "Cancer Type"] <- "Included as Cancer Type"
  
  return(GDSC_age)
}





############ anova_analysis_drugresponse_pan() ###################
# the 'anova_analysis_drugresponse_pan' function performs an ANOVA analysis to assess the effects of tissue/cancer type 
# and MSI status on drug response (LN_IC50). This ANOVA model helps determine whether different cancer types or MSI status 
# significantly influence drug sensitivity in cancer cell lines.
anova_analysis_drugresponse_pan <- function(GDSC_age) {
  # run ANOVA with tissue/cancer type and MSI status as factors
  anova_age_res <- aov(LN_IC50 ~ as.factor(tissue_or_cancer_chosen) + as.factor(msi_status), data = GDSC_age)
  anova_results_res <- summary(anova_age_res)
  
  GDSC_age$residuals <- residuals(anova_age_res)   # extract residuals
  
  # compute correlation between residuals and predicted age for each drug
  correlation_results_ANOVA <- GDSC_age %>%
    group_by(DRUG_NAME) %>%
    summarise(
      correlation = cor(residuals, age_prediction, use = "complete.obs"),
      p_value = cor.test(residuals, age_prediction, method = "pearson")$p.value
    ) %>%
    mutate(adj_p_value = p.adjust(p_value, method = "fdr")) %>%
    arrange(desc(correlation))
  
  # filter significant drugs
  significant_drugs <- correlation_results_ANOVA %>%
    filter(adj_p_value < 0.05)
  
  # count the number of significant drugs
  num_significant_drugs <- nrow(significant_drugs)
  print(paste("Number of drugs with adjusted p-value below", 0.05, ":", num_significant_drugs))
  
  return(list(
    anova_results = anova_results_res,
    residuals = GDSC_age$residuals,
    correlation_results = correlation_results_ANOVA,
    significant_drugs = significant_drugs
  ))
}






############ run_dsea_pan_cancer() ###################
# The 'run_dsea_pan_cancer' function performs a Drug Set Enrichment Analysis (DSEA) to assess whether there is still a signal
# in drug response data after ANOVA correction. This method ranks drugs based on the correlation between their
# predicted age and residual drug sensitivity (LN_IC50), then tests whether specific drug targets are significantly
# enriched among the top-ranked drugs.
run_dsea_pan_cancer <- function(correlation_results_ANOVA, GDSC_age_res) {
  set.seed(123)  
  
  # rank drugs based on correlation between predicted age and ANOVA residuals
  ranked_drugs <- correlation_results_ANOVA %>%
    arrange(desc(correlation))
  
  # create a named vector of rankings
  rankings_drugs <- setNames(ranked_drugs$correlation, ranked_drugs$DRUG_NAME)
  
  # define putative drug targets
  drug_targets <- GDSC_age_res %>%
    distinct(!!sym("DRUG_NAME"), .keep_all = TRUE) %>%
    separate_rows(!!sym("PUTATIVE_TARGET"), sep = ", ") %>%
    distinct(!!sym("DRUG_NAME"), !!sym("PUTATIVE_TARGET"))
  
  # create drug sets based on targets
  drug_sets <- split(drug_targets[["DRUG_NAME"]], drug_targets[["PUTATIVE_TARGET"]])
  
  # perform Drug Set Enrichment Analysis (DSEA)
  dsea_results <- fgsea(pathways = drug_sets, stats = rankings_drugs, minSize = 5, maxSize = 2000)
  
  # filter significant results
  significant_results_DSEA <- dsea_results %>%
    filter(padj < 0.05)
  
  # count significant results
  num_significant_targets <- nrow(significant_results_DSEA)
  print(paste("DSEA identified", num_significant_targets, "significantly enriched drug targets (padj < 0.05)."))
  print(significant_results_DSEA)
  
  significant_results_DSEA_df <- significant_results_DSEA %>%
    mutate(across(where(is.list), ~ sapply(., toString)))

  # save results 
  output_path <- "results/cell_lines_drug_sensitivity/dsea_pancancer/significant_enrichment_dsea_pancancer.xlsx"
  write_xlsx(significant_results_DSEA_df, output_path)
  
  return(significant_results_DSEA_df)
}





############ anova_drugresponse_cancer_specific() ###################
# the 'anova_drugresponse_cancer_specific' function performs an ANOVA analysis to assess the effect of MSI status 
# on drug response (LN_IC50) within each cancer type (when there is MSI). 
anova_drugresponse_cancer_specific <- function(GDSC_age_filtered) {
  GDSC_age_filtered$msi_status[GDSC_age_filtered$msi_status == ""] <- "MSS" # replace empty values with "MSS"
  
  unique_cancer_types <- unique(GDSC_age_filtered$cancer_type)
  
  anova_results_list <- list() # lists to store ANOVA and correlation results
  correlation_results_list <- list()
  
  cumulative_correlation_results <- data.frame() # data frame to store cumulative correlation results
  
  # loop through each cancer type
  for (cancer_type in unique_cancer_types) {
    cancer_data <- GDSC_age_filtered[GDSC_age_filtered$cancer_type == cancer_type, ] # filter data for the current cancer type
    
    if (length(unique(cancer_data$msi_status)) > 1) { # check if msi_status has at least two levels (excluding empty values)
      # perform ANOVA
      anova_model <- aov(LN_IC50 ~ as.factor(msi_status), data = cancer_data)
      
      cancer_data$residuals <- residuals(anova_model) # extract residuals
      
      # compute correlation between residuals and predicted age for each drug
      correlation_results <- cancer_data %>%
        group_by(DRUG_NAME) %>%
        summarise(
          correlation = if (n() > 2) cor(residuals, age_prediction, use = "complete.obs") else NA,
          p_value = if (n() > 2) cor.test(residuals, age_prediction, method = "pearson")$p.value else NA
        ) %>%
        mutate(adj_p_value = p.adjust(p_value, method = "fdr")) %>%
        arrange(desc(abs(correlation)))
      
      significant_drugs <- correlation_results %>% # filter significant results
        filter(!is.na(adj_p_value) & adj_p_value < 0.05)
      
      anova_results_list[[cancer_type]] <- list(anova_model = anova_model, significant_drugs = significant_drugs) # store results
    } else {
      # compute correlation between LN_IC50 and predicted age for each drug (without ANOVA)
      correlation_results <- cancer_data %>%
        group_by(DRUG_NAME) %>%
        summarise(
          correlation = if (n() > 2) cor(LN_IC50, age_prediction, use = "complete.obs") else NA,
          p_value = if (n() > 2) cor.test(LN_IC50, age_prediction, method = "pearson")$p.value else NA
        ) %>%
        mutate(adj_p_value = p.adjust(p_value, method = "fdr")) %>%
        arrange(desc(abs(correlation)))
      
      significant_drugs <- correlation_results %>% # filter significant results
        filter(!is.na(adj_p_value) & adj_p_value < 0.05)
      
      correlation_results_list[[cancer_type]] <- significant_drugs # store results
    }
    
    # append cumulative correlation results
    cumulative_correlation_results <- rbind(cumulative_correlation_results, 
                                            data.frame(cancer_type = cancer_type, correlation_results))
  }
  
  # print significant results for each cancer type
  for (cancer_type in names(anova_results_list)) {
    cat("Cancer Type:", cancer_type, "\n")
    print(anova_results_list[[cancer_type]]$significant_drugs)
  }
  
  for (cancer_type in names(correlation_results_list)) {
    cat("Cancer Type:", cancer_type, "\n")
    print(correlation_results_list[[cancer_type]])
  }
  
  return(list(
    anova_results = anova_results_list,
    correlation_results = correlation_results_list,
    cumulative_correlation_results = cumulative_correlation_results
  ))
}





############ summarize_significant_drugs() ###################
# the 'summarize_significant_drugs' function prints a summary of significant drugs per cancer type
# and saves them to an Excel file.
summarize_significant_drugs <- function(anova_results,
                                        output_file = "results/cell_lines_drug_sensitivity/significant_drugs_cancer_types/significant_drugs_anova_results.xlsx") {
  
  significant_summary <- "From the ANOVA analysis (cancer-type specific), the following significant drugs were identified:\n"
  
  # initialize an empty dataframe to store significant results
  significant_drugs_df <- data.frame(
    cancer_type = character(),
    DRUG_NAME = character(),
    correlation = numeric(),
    p_value = numeric(),
    adj_p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # loop through ANOVA results and collect significant drugs
  for (cancer_type in names(anova_results$anova_results)) {
    significant_drugs <- anova_results$anova_results[[cancer_type]]$significant_drugs
    if (!is.null(significant_drugs) && is.data.frame(significant_drugs) && nrow(significant_drugs) > 0) {
      drug_list <- unique(significant_drugs$DRUG_NAME)
      drug_list <- paste(drug_list, collapse = ", ")
      significant_summary <- paste(significant_summary, "-", cancer_type, ":", drug_list, "\n")
      
      # append to the dataframe
      significant_drugs_df <- bind_rows(significant_drugs_df, significant_drugs %>%
                                          mutate(cancer_type = cancer_type) %>%
                                          select(cancer_type, DRUG_NAME, correlation, p_value, adj_p_value))
    }
  }
  
  # loop through correlation results and collect significant drugs
  for (cancer_type in names(anova_results$correlation_results)) {
    significant_drugs <- anova_results$correlation_results[[cancer_type]]
    if (!is.null(significant_drugs) && is.data.frame(significant_drugs) && nrow(significant_drugs) > 0) {
      drug_list <- unique(significant_drugs$DRUG_NAME)
      drug_list <- paste(drug_list, collapse = ", ")
      significant_summary <- paste(significant_summary, "-", cancer_type, ":", drug_list, "\n")
      
      # append to the dataframe
      significant_drugs_df <- bind_rows(significant_drugs_df, significant_drugs %>%
                                          mutate(cancer_type = cancer_type) %>%
                                          select(cancer_type, DRUG_NAME, correlation, p_value, adj_p_value))
    }
  }
  
  # print summary
  cat(significant_summary)
  
  # save results to Excel only if significant drugs exist
  if (nrow(significant_drugs_df) > 0) {
    write_xlsx(significant_drugs_df, output_file)
    cat("\n Significant drug results saved to:", output_file, "\n")
  } else {
    cat("\n️ No significant drugs found. No file was saved.\n")
  }
}





############ extract_significant_drugs() ###################
# The 'extract_significant_drugs' function extract the significant drugs per cancer type
extract_significant_drugs <- function(anova_results) {
  significant_drugs_list <- list()
  
  # retrieve significant drugs from ANOVA results
  for (cancer_type in names(anova_results$anova_results)) {
    significant_drugs <- anova_results$anova_results[[cancer_type]]$significant_drugs
    if (!is.null(significant_drugs) && is.data.frame(significant_drugs) && nrow(significant_drugs) > 0) {
      significant_drugs_list[[cancer_type]] <- unique(significant_drugs$DRUG_NAME)
    }
  }
  
  # retrieve significant drugs from direct correlation results (for cancer types without ANOVA)
  for (cancer_type in names(anova_results$correlation_results)) {
    significant_drugs <- anova_results$correlation_results[[cancer_type]]
    if (!is.null(significant_drugs) && is.data.frame(significant_drugs) && nrow(significant_drugs) > 0) {
      significant_drugs_list[[cancer_type]] <- unique(significant_drugs$DRUG_NAME)
    }
  }
  
  print(significant_drugs_list)  # debugging step to check if all cancer types are present
  
  return(significant_drugs_list)
}





############ run_dsea_cancer_specific() ###################
# The 'run_dsea_cancer_specific' function performs Drug Set Enrichment Analysis (DSEA) for each cancer type.
# This method ranks drugs based on the correlation between predicted age and drug response (LN_IC50) 
# and tests whether specific drug targets are significantly enriched among the top-ranked drugs.
run_dsea_cancer_specific <- function(cumulative_correlation_results, GDSC_age_filtered) {
  set.seed(123)
  
  # filter correlations and compute mean per cancer type and drug
  correlation_results_filtered <- cumulative_correlation_results %>%
    filter(!is.na(correlation)) %>%
    group_by(cancer_type, DRUG_NAME) %>%
    summarise(correlation = mean(correlation, na.rm = TRUE), .groups = "drop")
  
  dsea_results_all <- list() # empty list to store results
  
  # extract acronyms for cancer types
  cancer_type_acronyms <- GDSC_age_filtered %>%
    select(cancer_type, acronyms_cancer_type) %>%
    distinct()
  
  # loop through each cancer type
  unique_cancer_types <- unique(correlation_results_filtered$cancer_type)
  
  for (cancer_type in unique_cancer_types) {
    # filter data for the current cancer type
    cancer_data <- correlation_results_filtered %>%
      filter(cancer_type == !!cancer_type)
    
    if (nrow(cancer_data) == 0) {
      print(paste("No data for cancer type:", cancer_type))
      next
    }
    
    # rank drugs for DSEA
    rankings <- setNames(cancer_data$correlation, cancer_data$DRUG_NAME)
    
    # define putative drug targets for the specific cancer type
    drug_targets <- GDSC_age_filtered %>%
      filter(cancer_type == !!cancer_type) %>%
      distinct(DRUG_NAME, PUTATIVE_TARGET) %>%
      separate_rows(PUTATIVE_TARGET, sep = ", ") %>%
      distinct(DRUG_NAME, PUTATIVE_TARGET)
    
    # create drug sets based on targets
    drug_sets <- split(drug_targets$DRUG_NAME, drug_targets$PUTATIVE_TARGET)
    
    # run Drug Set Enrichment Analysis (DSEA)
    dsea_results <- fgseaMultilevel(pathways = drug_sets, stats = rankings, minSize = 5, maxSize = 2000)
    
    # add cancer type and its acronym to the results
    acronym <- cancer_type_acronyms %>%
      filter(cancer_type == !!cancer_type) %>%
      pull(acronyms_cancer_type)
    
    dsea_results <- dsea_results %>%
      mutate(cancer_type = cancer_type, acronyms_cancer_type = acronym)
    
    dsea_results_all[[cancer_type]] <- dsea_results
  }
  
  dsea_combined_results <- bind_rows(dsea_results_all) # combine results into a single dataframe
  
  # significant results
  significant_results <- dsea_combined_results %>%
    filter(padj < 0.05)
  
  if (nrow(significant_results) > 0) {
    print("Significant drug targets found:")
    print(significant_results)
  } else {
    print("No significant targets found across all cancer types.")
  }
  
  # convert leading edge from list to string
  significant_results_df <- significant_results %>%
    mutate(across(where(is.list), ~ sapply(., toString)))
  
  output_path <- "results/cell_lines_drug_sensitivity/dsea_per_cancer_type/significant_enrichment_dsea_cancerspecific.xlsx"
  write_xlsx(significant_results_df, output_path)
  
  return(significant_results_df)
}





############ check_median_depletion() ###################
# the 'check_median_depletion' function is used to calculate and print median depletion scores for essential and non-essential genes
check_median_depletion <- function(scaled_DepMat, ess_genes, noness_genes) {
  
  # loop over each column (cell line) in the scaled_DepMat
  for (cell_line in colnames(scaled_DepMat)) {
    
    # extract the data for this cell line
    cell_line_data <- scaled_DepMat[, cell_line]
    
    # subset the essential and non-essential gene data for the current cell line
    ess_data <- cell_line_data[rownames(scaled_DepMat) %in% ess_genes]
    noness_data <- cell_line_data[rownames(scaled_DepMat) %in% noness_genes]
    
    # calculate the median depletion scores
    median_ess <- median(ess_data, na.rm = TRUE)
    median_noness <- median(noness_data, na.rm = TRUE)
    
    # print the results for this cell line
    cat("Cell Line:", cell_line, "\n")
    cat("Median Depletion for Essential Genes:", median_ess, "\n")
    cat("Median Depletion for Non-Essential Genes:", median_noness, "\n")
    cat("-------------------------------------------------\n")
  }
}


############ load_Rdata_files() ###################
# the 'load_Rdata' function is used to load all .RData files present in a specified directory
load_Rdata_files <- function(directory) {
  rdata_files <- list.files(directory, pattern = "\\.RData$", full.names = TRUE) # list of all .RData files in the specified directory
  
  for (file in rdata_files) {
    load(file, envir = .GlobalEnv)
    cat("Loaded:", file, "\n")  # print the file name to confirm it's loaded
  }
}





############ load_and_filter_genedep_pancancer() ###################
# the 'load_and_filter_genedep_pancancer' function processes gene dependency data
# from the DepMap dataset to identify essential genes across cell lines in a pan-cancer 
# context. The function applies several filtering steps to remove non-essential and 
# common fitness genes while retaining those significantly depleted in multiple cell lines.
#
# steps:
# 1. Load & Preprocess Data: 
#    - Reads the DepMap gene dependency matrix.
#    - Cleans gene names and transposes the matrix for easier manipulation.
#
# 2. Load Essential & Non-Essential Gene Lists: 
#    - Imports Achilles common essential genes and non-essential genes.
#    - Loads additional gene sets from the BAGEL dataset.
#
# 3. Scale Gene Dependency Scores: 
#    - Adjusts gene effect values relative to essential and non-essential gene distributions.
#    - Ensures the median depletion for essential genes is around -1 and for non-essential genes is around 0.
#
# 4. Identify & Remove Common Fitness Genes: 
#    - Uses the CoRe.FiPer algorithm to detect and exclude frequently depleted genes.
#    - Filters out additional essential genes related to key biological processes (e.g., DNA replication, ribosomal proteins).
#
# 5. Filter for Significant Dependencies: 
#    - Retains only genes that exhibit a depletion effect of -0.5 or lower in at least 'min_cell_lines' cell lines.
#
# 6. Final Processing & Alignment: 
#    - Matches the cell lines between DepMap and the provided cell line metadata (`cl_samples`).
#    - Ensures all genes and cell lines are aligned for downstream analysis.
#
# output:
# - a filtered and scaled gene dependency matrix (`gene_effect_filtered`).
#
# parameters:
# - depmap_file: path to the DepMap gene dependency dataset.
# - essential_genes_file: path to the Achilles Common Essential genes dataset.
# - nonessential_genes_file: path to the Achilles Never Essential genes dataset.
# - bagel_data_path: directory containing additional essential/non-essential gene lists from BAGEL.
# - cl_samples: data frame of cell line metadata for filtering.
# - min_cell_lines: minimum number of cell lines where a gene must show a depletion effect of -0.5 to be retained.
load_and_filter_genedep_pancancer <- function(depmap_file, essential_genes_file, nonessential_genes_file, 
                                              bagel_data_path, cl_samples, min_cell_lines) {
  
  library(tidyverse)
  library(CoRe)
  library(data.table)
  
  # ============================
  # LOAD GENE DEPENDENCY DATA
  # ============================
  
  # load gene effect data from DepMap
  DepMat <- fread(depmap_file)
  DepMat <- DepMat %>% column_to_rownames(var = "V1")
  
  # clean gene names 
  cleaned_gene_names <- gsub(" \\(.*\\)", "", colnames(DepMat))
  colnames(DepMat) <- cleaned_gene_names
  
  # transpose matrix to have genes as rows and cell lines as columns
  DepMat <- t(DepMat)
  
  # ============================
  # LOAD ESSENTIAL & NON-ESSENTIAL GENES
  # ============================
  
  # load Achilles 2024 essential and non essential genes
  common_essential_genes <- read_csv(essential_genes_file)
  non_essential_genes <- read_csv(nonessential_genes_file)
  
  # clean names and create vectors
  ess_genes <- gsub(" \\(.*\\)", "", common_essential_genes$Gene)
  noness_genes <- gsub(" \\(.*\\)", "", non_essential_genes$Gene)
  
  # ============================
  # SCALE DEPENDENCY MATRIX
  # ============================
  
  # scale using essential and non-essential genes 
  scaled_DepMat <- CoRe.scale_to_essentials(DepMat,
                                            ess_genes,
                                            noness_genes)
  
  # validate scaling (should return 0 for non-essentials and -1 for essentials)
  check_median_depletion(scaled_DepMat, ess_genes, noness_genes)
  
  # identify core fitness genes using FiPer
  scaled_DepMat_Fiper <- CoRe.FiPer(scaled_DepMat,
                                    display=TRUE,
                                    percentile=0.9,
                                    method='AUC')
  # extract the found core fitness genes
  CoRefitness_genes <- scaled_DepMat_Fiper$cfgenes
  
  # ============================
  # LOAD ADDITIONAL FILTERING GENES
  # ============================
  
  # load additional essential genes (BAGEL genes, CURATED BAGEL genes, ess. genes DNA REPLICATION,
  # ess. genes HISTONES, ess. genes KEGG RNA polymerase, ess. genes PROTEASOME,
  # ess. genes RIBOSOMAL PROTEINS, ess. genes SPLICEOSOME)
  load_Rdata_files(bagel_data_path)
  
  # combine all genes to filter
  genes_to_filter <- unique(c(CoRefitness_genes, 
                              BAGEL_essential, 
                              BAGEL_nonEssential,
                              EssGenes.DNA_REPLICATION_cons,
                              EssGenes.HISTONES, 
                              EssGenes.KEGG_rna_polymerase,
                              EssGenes.PROTEASOME_cons, 
                              EssGenes.ribosomalProteins,
                              EssGenes.SPLICEOSOME_cons, 
                              ess_genes,
                              noness_genes))
  
  # remove filtered genes
  filtered_genes <- rownames(scaled_DepMat)[rownames(scaled_DepMat) %in% genes_to_filter]
  gene_effect <- scaled_DepMat[!(rownames(scaled_DepMat) %in% filtered_genes), ]
  
  print(dim(gene_effect))  # check gene count
  write.csv(gene_effect, "data/depmap/filtered_scaled_DepMat_2024.csv", row.names = TRUE)
  
  # ============================
  # MATCHING CELL LINES
  # ============================
  
  # load processed dependency matrix
  matching_cell_lines <- intersect(colnames(gene_effect), cl_samples$BROAD_ID)
  print(paste("Number of matching cell lines:", length(matching_cell_lines)))
  
  # filter cell lines shared by Sanger and DepMap
  cl_samples_gene_filtered <- cl_samples[cl_samples$BROAD_ID %in% matching_cell_lines, ]
  gene_effect_filtered <- gene_effect[, colnames(gene_effect) %in% matching_cell_lines]
  
  # remove genes with NA values
  gene_effect_filtered <- gene_effect_filtered[complete.cases(gene_effect_filtered), ]
  print(dim(gene_effect_filtered))  # check final dimensions
  
  # check the order
  cl_samples_gene_filtered <- cl_samples_gene_filtered[match(colnames(gene_effect_filtered), cl_samples_gene_filtered$BROAD_ID), ]
  if (!all(colnames(gene_effect_filtered) == cl_samples_gene_filtered$BROAD_ID)) {
    stop("Error: Cell lines are not aligned.")
  } else {
    print("Cell lines are correctly aligned between gene_effect and cl_samples.")
  }
  
  # ============================
  # FILTER SIGNIFICANT GENES
  # ============================
  
  # identify genes with depletion effect <= -0.5 in at least 'min_cell_lines' cell lines
  gene_counts <- apply(gene_effect_filtered, 1, function(x) sum(x <= -0.5, na.rm = TRUE))
  essential_genes <- names(gene_counts[gene_counts >= min_cell_lines])
  
  # keep only significant genes
  gene_effect_filtered <- gene_effect_filtered[rownames(gene_effect_filtered) %in% essential_genes, ]
  
  # transpose matrix for easier downstream analysis
  gene_effect_filtered <- t(gene_effect_filtered)
  print(dim(gene_effect_filtered))  # final dimensions
  
  return(list(gene_effect_filtered = gene_effect_filtered, 
               cl_samples_gene_filtered = cl_samples_gene_filtered))
}





############ generate_merged_gene_dataset() ###################
# The 'generate_merged_gene_dataset' function filters a dataset based on cancer type or tissue availability, selecting 
# categories with a specified minimum number of cell lines and classifying each entry by its cancer type or tissue.

generate_merged_gene_dataset <- function(gene_effect_filtered, cl_samples, threshold) {   
  # merge gene dependency data with cell line metadata using BROAD_ID
  depmap_combined <- left_join(
    cl_samples[, c("BROAD_ID", "age_prediction", "age_group_median", "age_at_sampling",
                   "tissue", "msi_status", "cancer_type", "acronyms_cancer_type", "CAccession")], 
    as.data.frame(gene_effect_filtered) %>% rownames_to_column("BROAD_ID"),
    by = "BROAD_ID"
  )
  
  # remove rows with missing values in key metadata fields
  depmap_combined <- depmap_combined[!is.na(depmap_combined$age_prediction) & 
                                       !is.na(depmap_combined$cancer_type) & 
                                       !is.na(depmap_combined$msi_status), ]
  
  print(dim(depmap_combined)) # check dimensions after merging
  
  # apply tissue or cancer filtering function
  depmap_combined <- filter_tissue_or_cancer(depmap_combined, "cancer_type", "tissue", "BROAD_ID", threshold)
  
  # update missing values in classification column
  depmap_combined$tissue_or_cancer_chosen[is.na(depmap_combined$tissue_or_cancer_chosen)] <- "others"
  
  # create 'status' column for summary plots
  depmap_combined$status <- ifelse(depmap_combined$tissue_or_cancer_chosen == "others", "Others", 
                                   ifelse(depmap_combined$tissue_or_cancer_chosen %in% depmap_combined$cancer_type, "Cancer Type", "Tissue Type"))
  
  # define a column for cancer-specific analysis
  depmap_combined$cancer_specific <- "Excluded"
  
  # mark samples as "Included as Cancer Type" if they belong to a sufficiently large cancer type
  depmap_combined$cancer_specific[depmap_combined$status == "Cancer Type"] <- "Included as Cancer Type"
  
  # reorder columns for clarity
  col_order <- colnames(depmap_combined)
  new_order <- c(col_order[1:9], "status", "tissue_or_cancer_chosen", col_order[10:ncol(depmap_combined)])
  depmap_combined <- depmap_combined[, new_order]
  
  # select only numeric columns (genes) for pivot_longer
  cols_to_exclude <- c("BROAD_ID", "age_prediction", "age_group_median", "age_at_sampling",
                       "tissue", "msi_status", "cancer_type", "acronyms_cancer_type", "CAccession",
                       "cancer_specific", "status", "tissue_or_cancer_chosen")
  
  numeric_cols <- depmap_combined %>% select(-all_of(cols_to_exclude)) %>% select(where(is.numeric)) %>% names()
  
  # pivot the data into long format
  gene_dep_age <- pivot_longer(
    depmap_combined,
    cols = all_of(numeric_cols),  # Use only numeric columns
    names_to = "Gene",  # Column for gene names
    values_to = "gene_effect"  # Column for gene effect values
  )
  
  # remove duplicates 
  gene_dep_age$status.1 <- NULL
  gene_dep_age$tissue_or_cancer_chosen.1 <- NULL
  
  print(dim(gene_dep_age)) # check final dimensions
  
  return(list(gene_dep_age = gene_dep_age, 
              depmap_combined = depmap_combined))
}





############ anova_analysis_genedep_pan() ###################
# The 'anova_analysis_genedep_pan' function performs an ANOVA analysis to assess 
# the effects of tissue/cancer type and MSI status on gene dependency (gene effect). 
# It then computes the correlation between the residuals and predicted age.

anova_analysis_genedep_pan <- function(gene_dep_age) {   
  # run ANOVA with tissue/cancer type and MSI status as factors
  anova_age_res <- aov(gene_effect ~ as.factor(tissue_or_cancer_chosen) + as.factor(msi_status), data = gene_dep_age)   
  anova_results_res <- summary(anova_age_res)  # store ANOVA results
  
  # extract residuals from the model
  gene_dep_age$residuals <- residuals(anova_age_res)  
  
  # compute correlation between residuals and predicted age for each gene
  correlation_results_ANOVA <- gene_dep_age %>% 
    group_by(Gene) %>% 
    summarise(
      correlation = cor(residuals, age_prediction, use = "complete.obs"),
      p_value = cor.test(residuals, age_prediction, method = "pearson")$p.value
    ) %>% 
    mutate(adj_p_value = p.adjust(p_value, method = "fdr")) %>% 
    arrange(desc(correlation))
  
  # filter significant genes
  significant_genes <- correlation_results_ANOVA %>% 
    filter(adj_p_value < 0.05)  
  
  # count the number of significant genes
  num_significant_genes <- nrow(significant_genes)
  print(paste("Number of genes with adjusted p-value below", 0.05, ":", num_significant_genes))
  
  return(list(
    anova_results = anova_results_res,
    residuals = gene_dep_age$residuals,
    correlation_results = correlation_results_ANOVA,
    significant_genes = significant_genes
  ))
}







############ run_gsea_pan_cancer() ###################
# the 'run_gsea_pan_cancer' function performs a Gene Set Enrichment Analysis (GSEA) 
# to assess whether there is still a biological signal in gene dependency data after ANOVA correction.
# this method ranks genes based on the correlation between their predicted age and ANOVA residuals,
# then tests whether specific biological pathways (GO BP, GO MF, GO CC, KEGG) are significantly 
# enriched among the top-ranked genes.

run_gsea_pan_cancer <- function(correlation_results_ANOVA) {   
  set.seed(123)  
  
  # rank genes based on correlation between predicted age and ANOVA residuals
  ranked_genes <- correlation_results_ANOVA %>% 
    arrange(desc(correlation))  # sort in descending order
  
  # create a named vector of rankings
  rankings <- setNames(ranked_genes$correlation, ranked_genes$Gene)
  
  # convert gene symbols to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db, 
                       keys = names(rankings), 
                       column = "ENTREZID", 
                       keytype = "SYMBOL", 
                       multiVals = "first")
  
  entrez_ids <- na.omit(entrez_ids)
  
  # assign Entrez IDs to rankings
  names(rankings) <- entrez_ids[names(rankings)]
  
  # ============================
  # GSEA ENRICHMENT ANALYSIS
  # ============================
  
  # GO Biological Process (BP)
  gseabp <- gseGO(
    geneList = rankings, 
    OrgDb = org.Hs.eg.db, 
    keyType = "ENTREZID", 
    ont = "BP", 
    minGSSize = 5, 
    maxGSSize = 2000, 
    pvalueCutoff = 0.05
  )
  
  # GO Molecular Function (MF)
  gseamf <- gseGO(
    geneList = rankings, 
    OrgDb = org.Hs.eg.db, 
    keyType = "ENTREZID", 
    ont = "MF", 
    minGSSize = 5, 
    maxGSSize = 2000, 
    pvalueCutoff = 0.05
  )
  
  # GO Cellular Component (CC)
  gseacc <- gseGO(
    geneList = rankings, 
    OrgDb = org.Hs.eg.db, 
    keyType = "ENTREZID", 
    ont = "CC", 
    minGSSize = 5, 
    maxGSSize = 2000, 
    pvalueCutoff = 0.05
  )
  
  # KEGG Pathways
  gseakegg <- gseKEGG(
    geneList = rankings, 
    organism = "hsa", 
    minGSSize = 5, 
    maxGSSize = 2000, 
    pvalueCutoff = 0.05
  )
  
  # ============================
  # EXTRACT SIGNIFICANT RESULTS
  # ============================
  significant_gseabp <- gseabp@result %>% filter(p.adjust < 0.05)
  significant_gseamf <- gseamf@result %>% filter(p.adjust < 0.05)
  significant_gseacc <- gseacc@result %>% filter(p.adjust < 0.05)
  significant_gseakegg <- gseakegg@result %>% filter(p.adjust < 0.05)
  
  # number of significant pathways found
  print(paste("Significant GO BP terms:", nrow(significant_gseabp)))
  print(paste("Significant GO MF terms:", nrow(significant_gseamf)))
  print(paste("Significant GO CC terms:", nrow(significant_gseacc)))
  print(paste("Significant KEGG pathways:", nrow(significant_gseakegg)))
  
  # ============================
  # SAVE SIGNIFICANT RESULTS
  # ============================
  
  # remove "core_enrichment" column before saving
  significant_gseabp_df <- significant_gseabp[, !colnames(significant_gseabp) %in% "core_enrichment"]
  significant_gseamf_df <- significant_gseamf[, !colnames(significant_gseamf) %in% "core_enrichment"]
  significant_gseacc_df <- significant_gseacc[, !colnames(significant_gseacc) %in% "core_enrichment"]
  significant_gseakegg_df <- significant_gseakegg[, !colnames(significant_gseakegg) %in% "core_enrichment"]
  
  write_xlsx(significant_gseabp_df, "results/cell_lines_gene_dependencies/gsea_pancancer/GOBP_enrichments.xlsx")
  write_xlsx(significant_gseamf_df, "results/cell_lines_gene_dependencies/gsea_pancancer/GOMF_enrichments.xlsx")
  write_xlsx(significant_gseacc_df, "results/cell_lines_gene_dependencies/gsea_pancancer/GOCC_enrichments.xlsx")
  write_xlsx(significant_gseakegg_df, "results/cell_lines_gene_dependencies/gsea_pancancer/GOKEGG_enrichments.xlsx")
  
  return(list(
    GO_BP = significant_gseabp_df,
    GO_MF = significant_gseamf_df,
    GO_CC = significant_gseacc_df,
    KEGG = significant_gseakegg_df,
    rankings = rankings
  ))
}



############ load_and_filter_genedep_cancer_type() ###################
# the 'load_and_filter_genedep_cancer_type' function processes gene dependency data
# from the DepMap dataset to identify essential genes across cell lines in a cancer specific
# context. The function applies several filtering steps to remove non-essential and 
# common fitness genes while retaining those significantly depleted in multiple cell lines.
load_and_filter_genedep_cancer_type <- function(depmap_file, essential_genes_file, nonessential_genes_file, 
                                                bagel_data_path, cl_samples, min_cell_lines = 2, thresh_number = 10) {
  
  library(tidyverse)
  library(CoRe)
  library(data.table)
  
  # ============================
  # LOAD GENE DEPENDENCY DATA
  # ============================
  
  DepMat <- fread(depmap_file) %>% column_to_rownames(var = "V1")
  colnames(DepMat) <- gsub(" \\(.*\\)", "", colnames(DepMat))  # Clean gene names
  DepMat <- t(DepMat)  # Transpose for easier manipulation
  
  # ============================
  # LOAD ESSENTIAL & NON-ESSENTIAL GENES
  # ============================
  
  common_essential_genes <- read_csv(essential_genes_file)
  non_essential_genes <- read_csv(nonessential_genes_file)
  
  ess_genes <- gsub(" \\(.*\\)", "", common_essential_genes$Gene)
  noness_genes <- gsub(" \\(.*\\)", "", non_essential_genes$Gene)
  
  # ============================
  # MATCHING CELL LINES
  # ============================
  
  matching_cell_lines <- intersect(colnames(DepMat), cl_samples$BROAD_ID)
  cl_samples_filtered <- cl_samples[cl_samples$BROAD_ID %in% matching_cell_lines, ]
  DepMat_filtered <- DepMat[, colnames(DepMat) %in% matching_cell_lines]
  
  cl_samples_filtered <- cl_samples_filtered[match(colnames(DepMat_filtered), cl_samples_filtered$BROAD_ID), ]
  if (!all(colnames(DepMat_filtered) == cl_samples_filtered$BROAD_ID)) {
    stop("Error: Cell lines are not aligned.")
  }
  
  # ============================
  # FILTER CANCER TYPES
  # ============================
  
  cancer_types_to_keep <- tapply(cl_samples_filtered$BROAD_ID, cl_samples_filtered$cancer_type, function(i) length(unique(i))) >= thresh_number
  valid_cancer_types <- names(cancer_types_to_keep)[cancer_types_to_keep]
  
  scaled_DepMat_list <- list()
  
  for (cancer_type in valid_cancer_types) {
    
    cancer_data <- cl_samples_filtered[cl_samples_filtered$cancer_type == cancer_type, ]
    DepMat_cancer <- DepMat_filtered[, colnames(DepMat_filtered) %in% cancer_data$BROAD_ID]
    
    # Load BAGEL and curated essential genes
    load_Rdata_files(bagel_data_path)
    genes_to_filter <- unique(c(ess_genes, noness_genes, BAGEL_essential, BAGEL_nonEssential,
                                EssGenes.DNA_REPLICATION_cons, EssGenes.HISTONES,
                                EssGenes.KEGG_rna_polymerase, EssGenes.PROTEASOME_cons,
                                EssGenes.ribosomalProteins, EssGenes.SPLICEOSOME_cons))
    
    # Scale using CoRe
    scaled_DepMat_cancer <- CoRe.scale_to_essentials(DepMat_cancer, ess_genes, noness_genes)
    CoRefitness_genes <- CoRe.FiPer(scaled_DepMat_cancer, display = TRUE, percentile = 0.9, method = 'AUC')$cfgenes
    
    # Remove genes in genes_to_filter
    genes_to_filter <- unique(c(genes_to_filter, CoRefitness_genes))
    filtered_scaled_DepMat_cancer <- scaled_DepMat_cancer[!(rownames(scaled_DepMat_cancer) %in% genes_to_filter), ]
    filtered_scaled_DepMat_cancer <- filtered_scaled_DepMat_cancer[complete.cases(filtered_scaled_DepMat_cancer), ]
    
    # Keep genes with a median depletion <= -0.5 in at least `min_cell_lines`
    gene_counts <- apply(filtered_scaled_DepMat_cancer, 1, function(x) sum(x <= -0.5, na.rm = TRUE))
    essential_genes <- names(gene_counts[gene_counts >= min_cell_lines])
    filtered_scaled_DepMat_cancer <- filtered_scaled_DepMat_cancer[rownames(filtered_scaled_DepMat_cancer) %in% essential_genes, ]
    
    scaled_DepMat_list[[cancer_type]] <- filtered_scaled_DepMat_cancer
  }
  
  # save scaled and filtered matrices for each cancer type
  for (cancer_type in names(scaled_DepMat_list)) {
    write.csv(scaled_DepMat_list[[cancer_type]], 
              file = paste0("data/depmap/cancer_specific/scaled_filtered_DepMat_", cancer_type, "_2024.csv"), 
              row.names = TRUE)
  }
    
    return(scaled_DepMat_list)
  
}
  




############ anova_genedep_cancer_specific() ###################
# the 'anova_genedep_cancer_specific' function performs an ANOVA analysis to assess the effect of MSI status
# on gene dependency scores within each cancer type (when MSI is present).
# if MSI status is not available, a direct correlation analysis between gene effect and predicted age is performed.
anova_genedep_cancer_specific <- function(scaled_DepMat_list, cl_samples_filtered) {
  
  # replace empty MSI values with "MSS" for consistency
  cl_samples_filtered$msi_status[cl_samples_filtered$msi_status == ""] <- "MSS"
  
  # lists to store ANOVA and correlation results
  anova_results_list <- list()
  correlation_results_list <- list()
  cumulative_correlation_results <- data.frame()
  
  # loop through each cancer type
  for (cancer_type in names(scaled_DepMat_list)) {
    
    # retrieve the gene dependency matrix for the specific cancer type
    gene_effect_filtered <- scaled_DepMat_list[[cancer_type]]
    
    # extract metadata for the selected cancer type
    cancer_data <- cl_samples_filtered[cl_samples_filtered$cancer_type == cancer_type, ]
    
    # check if gene dependency matrix and metadata are aligned
    if (!all(colnames(gene_effect_filtered) == cancer_data$BROAD_ID)) {
      stop(paste("Error: cell lines are not aligned for", cancer_type))
    }
    
    # transpose the matrix for analysis
    gene_effect_filtered <- t(gene_effect_filtered)
    
    # merge gene dependency data with metadata
    merged_data <- cbind(cancer_data[, c("BROAD_ID", "tissue", "cancer_type","acronyms_cancer_type", "msi_status", 
                                         "age_prediction", "age_group_median")], 
                         gene_effect_filtered)
    
    
    # convert data to long format (each row corresponds to a gene and a cell line)
    gene_dep_age <- pivot_longer(merged_data, 
                                 cols = 8:ncol(merged_data),  
                                 names_to = "Gene", 
                                 values_to = "gene_effect")
    
    # check if MSI status has more than one level, enabling ANOVA
    if (length(unique(gene_dep_age$msi_status)) > 1) {
      
      # perform ANOVA to remove variability due to MSI status
      anova_age_res <- aov(gene_effect ~ as.factor(msi_status), data = gene_dep_age)
      
      # extract residuals from the ANOVA model
      residuals_anova <- residuals(anova_age_res)
      
      # assign residuals only to valid rows (non-missing MSI status and gene effect)
      valid_rows <- !is.na(gene_dep_age$msi_status) & !is.na(gene_dep_age$gene_effect)
      gene_dep_age$residuals[valid_rows] <- residuals_anova
      
      # filter out incomplete rows
      gene_dep_age <- gene_dep_age[complete.cases(gene_dep_age[, c("residuals", "age_prediction")]), ]
      
      # compute correlation between residuals and predicted age for each gene
      correlation_results <- gene_dep_age %>%
        group_by(Gene) %>%
        summarise(
          correlation = if (n() > 2) cor(residuals, age_prediction, use = "complete.obs") else NA,
          p_value = if (n() > 2) cor.test(residuals, age_prediction, method = "pearson")$p.value else NA
        ) %>%
        mutate(adj_p_value = p.adjust(p_value, method = "fdr")) %>%
        arrange(desc(abs(correlation)))
      
      # filter significant genes
      significant_genes <- correlation_results %>%
        filter(!is.na(adj_p_value) & adj_p_value < 0.05) %>%
        arrange(desc(abs(correlation)))
      
      # store ANOVA results and significant genes
      anova_results_list[[cancer_type]] <- list(anova_model = anova_age_res, significant_genes = significant_genes)
      
    } else {
      
      # if MSI is not applicable, perform direct correlation between gene effect and predicted age
      gene_dep_age <- gene_dep_age[complete.cases(gene_dep_age[, c("gene_effect", "age_prediction")]), ]
      
      correlation_results <- gene_dep_age %>%
        group_by(Gene) %>%
        summarise(
          correlation = if (n() > 2) cor(gene_effect, age_prediction, use = "complete.obs") else NA,
          p_value = if (n() > 2) cor.test(gene_effect, age_prediction, method = "pearson")$p.value else NA
        ) %>%
        mutate(adj_p_value = p.adjust(p_value, method = "fdr")) %>%
        arrange(desc(abs(correlation)))
      
      # filter significant genes
      significant_genes <- correlation_results %>%
        filter(!is.na(adj_p_value) & adj_p_value < 0.05) %>%
        arrange(desc(abs(correlation)))
      
      # store correlation results
      correlation_results_list[[cancer_type]] <- significant_genes
    }
    
    # add cumulative correlation results
    cumulative_correlation_results <- rbind(cumulative_correlation_results, 
                                            data.frame(cancer_type = cancer_type,
                                                       acronyms_cancer_type = unique(cancer_data$acronyms_cancer_type),
                                                       correlation_results))
  }
  
  # print significant results for each cancer type
  for (cancer_type in names(anova_results_list)) {
    cat("Cancer Type:", cancer_type, "\n")
    print(anova_results_list[[cancer_type]]$significant_genes)
  }
  
  for (cancer_type in names(correlation_results_list)) {
    cat("Cancer Type:", cancer_type, "\n")
    print(correlation_results_list[[cancer_type]])
  }
  
  # return results
  return(list(
    anova_results = anova_results_list,
    correlation_results = correlation_results_list,
    cumulative_correlation_results = cumulative_correlation_results
  ))
}





############ summarize_significant_genes() ###################
# the 'summarize_significant_genes' function prints a summary of significant genes per cancer type
# and saves them to an Excel file.
summarize_significant_genes <- function(anova_results, 
                                        output_file = "results/cell_lines_gene_dependencies/significant_genes_cancer_types/significant_genes_anova_results.xlsx") {
  
  significant_summary <- "From the ANOVA analysis (cancer-type specific), the following significant genes were identified:\n"
  
  # initialize an empty dataframe to store significant results
  significant_genes_df <- data.frame(
    cancer_type = character(),
    Gene = character(),
    correlation = numeric(),
    p_value = numeric(),
    adj_p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # loop through ANOVA results and collect significant genes
  for (cancer_type in names(anova_results$anova_results)) {
    significant_genes <- anova_results$anova_results[[cancer_type]]$significant_genes
    if (!is.null(significant_genes) && is.data.frame(significant_genes) && nrow(significant_genes) > 0) {
      gene_list <- unique(significant_genes$Gene)
      gene_list <- paste(gene_list, collapse = ", ")
      significant_summary <- paste(significant_summary, "-", cancer_type, ":", gene_list, "\n")
      
      # append to the dataframe
      significant_genes_df <- bind_rows(significant_genes_df, significant_genes %>%
                                          mutate(cancer_type = cancer_type) %>%
                                          select(cancer_type, Gene, correlation, p_value, adj_p_value))
    }
  }
  
  # loop through correlation results and collect significant genes
  for (cancer_type in names(anova_results$correlation_results)) {
    significant_genes <- anova_results$correlation_results[[cancer_type]]
    if (!is.null(significant_genes) && is.data.frame(significant_genes) && nrow(significant_genes) > 0) {
      gene_list <- unique(significant_genes$Gene)
      gene_list <- paste(gene_list, collapse = ", ")
      significant_summary <- paste(significant_summary, "-", cancer_type, ":", gene_list, "\n")
      
      # append to the dataframe
      significant_genes_df <- bind_rows(significant_genes_df, significant_genes %>%
                                          mutate(cancer_type = cancer_type) %>%
                                          select(cancer_type, Gene, correlation, p_value, adj_p_value))
    }
  }
  
  # print summary
  cat(significant_summary)
  
  # save results to Excel only if significant genes exist
  if (nrow(significant_genes_df) > 0) {
    write_xlsx(significant_genes_df, output_file)
    cat("\n Significant gene results saved to:", output_file, "\n")
  } else {
    cat("\n⚠ No significant genes found. No file was saved.\n")
  }
}





############ extract_significant_genes() ###################
# the 'extract_significant_genes' function extract the significant genes per cancer type
extract_significant_genes <- function(anova_results) {
  
  significant_genes_list <- list() # initialize list
  
  # retrieve significant genes from ANOVA results
  for (cancer_type in names(anova_results$anova_results)) {
    significant_genes <- anova_results$anova_results[[cancer_type]]$significant_genes
    
    if (!is.null(significant_genes) && is.data.frame(significant_genes) && nrow(significant_genes) > 0) {
      significant_genes_list[[cancer_type]] <- unique(significant_genes$Gene)
    }
  }
  
  # retrieve significant genes from direct correlation results (for cancer types without ANOVA)
  for (cancer_type in names(anova_results$correlation_results)) {
    significant_genes <- anova_results$correlation_results[[cancer_type]]
    
    if (!is.null(significant_genes) && is.data.frame(significant_genes) && nrow(significant_genes) > 0) {
      significant_genes_list[[cancer_type]] <- unique(significant_genes$Gene)
    }
  }
  
  print(significant_genes_list)  # debugging step to check extracted genes
  
  return(significant_genes_list)
}





############ run_gsea_cancer_specific() ###################
# The 'run_gsea_cancer_specific' function performs Gene Set Enrichment Analysis (GSEA) for each cancer type.
# This method ranks grnrd based on the correlation between predicted age and gene dependency (gene effect) 
# and tests whether specific pathways are significantly enriched among the top-ranked genes.
run_gsea_cancer_specific <- function(cumulative_correlation_results, 
                                     output_dir = "results/cell_line_gene_dependencies/gsea_per_cancer_type") {
  
  set.seed(123)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE) # create output directory if it does not exist
  
  # filter correlations and compute mean per cancer type and gene
  correlation_results_filtered <- cumulative_correlation_results %>%
    filter(!is.na(correlation)) %>%
    group_by(cancer_type, Gene) %>%
    summarise(correlation = mean(correlation, na.rm = TRUE), .groups = "drop")
  
  gsea_results_all <- list() # initialize empty list
  
  unique_cancer_types <- unique(correlation_results_filtered$cancer_type) # get unique cancer types
  
  for (cancer_type in unique_cancer_types) {
    
    print(paste("Running GSEA for", cancer_type, "..."))
    
    # filter data for the current cancer type
    cancer_data <- correlation_results_filtered %>%
      filter(cancer_type == !!cancer_type)
    
    if (nrow(cancer_data) == 0) {
      print(paste("No data for cancer type:", cancer_type))
      next
    }
    
    # create rankings vector
    rankings <- setNames(cancer_data$correlation, cancer_data$Gene)
    
    # convert gene symbols to ENTREZ IDs
    entrez_ids <- bitr(names(rankings), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    if (nrow(entrez_ids) < length(rankings)) {
      warning(paste0(round((1 - nrow(entrez_ids) / length(rankings)) * 100, 2),
                     "% of genes could not be mapped to ENTREZ IDs."))
    }
    
    rankings <- rankings[names(rankings) %in% entrez_ids$SYMBOL]
    names(rankings) <- entrez_ids$ENTREZID[match(names(rankings), entrez_ids$SYMBOL)]
    
    # sort rankings in decreasing order
    rankings <- rankings[order(-rankings)]
    
    # run GSEA for different gene sets
    gseabp <- gseGO(geneList = rankings, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                    minGSSize = 5, maxGSSize = 2000, pvalueCutoff = 0.05)
    
    gseamf <- gseGO(geneList = rankings, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF",
                    minGSSize = 5, maxGSSize = 2000, pvalueCutoff = 0.05)
    
    gseacc <- gseGO(geneList = rankings, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC",
                    minGSSize = 5, maxGSSize = 2000, pvalueCutoff = 0.05)
    
    gseakegg <- gseKEGG(geneList = rankings, organism = "hsa", minGSSize = 5, maxGSSize = 2000, pvalueCutoff = 0.05)
    
    # store results
    gsea_results_all[[cancer_type]] <- list(
      BP = as.data.frame(gseabp@result) %>% mutate(cancer_type = cancer_type),
      MF = as.data.frame(gseamf@result) %>% mutate(cancer_type = cancer_type),
      CC = as.data.frame(gseacc@result) %>% mutate(cancer_type = cancer_type),
      KEGG = as.data.frame(gseakegg@result) %>% mutate(cancer_type = cancer_type)
    )
    
    # save individual results per cancer type
    write.csv(gsea_results_all[[cancer_type]]$BP, file = file.path(output_dir, paste0(cancer_type, "_BP_results.csv")), row.names = FALSE)
    write.csv(gsea_results_all[[cancer_type]]$MF, file = file.path(output_dir, paste0(cancer_type, "_MF_results.csv")), row.names = FALSE)
    write.csv(gsea_results_all[[cancer_type]]$CC, file = file.path(output_dir, paste0(cancer_type, "_CC_results.csv")), row.names = FALSE)
    write.csv(gsea_results_all[[cancer_type]]$KEGG, file = file.path(output_dir, paste0(cancer_type, "_KEGG_results.csv")), row.names = FALSE)
    
  }
  
  # combine results across all cancer types
  combined_results_bp <- do.call(rbind, lapply(gsea_results_all, `[[`, "BP"))
  combined_results_mf <- do.call(rbind, lapply(gsea_results_all, `[[`, "MF"))
  combined_results_cc <- do.call(rbind, lapply(gsea_results_all, `[[`, "CC"))
  combined_results_kegg <- do.call(rbind, lapply(gsea_results_all, `[[`, "KEGG"))
  
  # save combined results
  write.csv(combined_results_bp, file = file.path(output_dir, "01_combined_BP_results.csv"), row.names = FALSE)
  write.csv(combined_results_mf, file = file.path(output_dir, "02_combined_MF_results.csv"), row.names = FALSE)
  write.csv(combined_results_cc, file = file.path(output_dir, "03_combined_CC_results.csv"), row.names = FALSE)
  write.csv(combined_results_kegg, file = file.path(output_dir, "04_combined_KEGG_results.csv"), row.names = FALSE)
  
  # filter significant results (FDR < 0.05)
  significant_results_bp <- combined_results_bp %>% filter(p.adjust < 0.05)
  significant_results_mf <- combined_results_mf %>% filter(p.adjust < 0.05)
  significant_results_cc <- combined_results_cc %>% filter(p.adjust < 0.05)
  significant_results_kegg <- combined_results_kegg %>% filter(p.adjust < 0.05)
  
  # print significant results
  print("Significant BP Results")
  print(head(significant_results_bp))
  
  print("Significant MF Results")
  print(head(significant_results_mf))
  
  print("Significant CC Results")
  print(head(significant_results_cc))
  
  print("Significant KEGG Results")
  print(head(significant_results_kegg))
  
  # save significant results (excluding 'core_enrichment' column)
  significant_results_bp_df <- significant_results_bp[, !colnames(significant_results_bp) %in% "core_enrichment"]
  significant_results_mf_df <- significant_results_mf[, !colnames(significant_results_mf) %in% "core_enrichment"]
  significant_results_cc_df <- significant_results_cc[, !colnames(significant_results_cc) %in% "core_enrichment"]
  significant_results_kegg_df <- significant_results_kegg[, !colnames(significant_results_kegg) %in% "core_enrichment"]
  
  write_xlsx(significant_results_bp_df, file.path(output_dir, "combined_GOBP_enrichments.xlsx"))
  write_xlsx(significant_results_mf_df, file.path(output_dir, "combined_GOMF_enrichments.xlsx"))
  write_xlsx(significant_results_cc_df, file.path(output_dir, "combined_GOCC_enrichments.xlsx"))
  write_xlsx(significant_results_kegg_df, file.path(output_dir, "combined_KEGG_enrichments.xlsx"))
  
  return(list(
    BP = significant_results_bp,
    MF = significant_results_mf,
    CC = significant_results_cc,
    KEGG = significant_results_kegg
  ))
}



























############ anova_exp_settings_analysis() ###################
# The 'anova_exp_settings_analysis' function performs an ANOVA model
# to adjust for cancer type/tissue and MSI status on a selected experimental response variable.
# Then, it calculates the correlation between the residuals and a selected age-related variable.
# arguments:
# - cl_samples: a dataframe containing the required columns
# - response_variable: the name of the experimental response variable (e.g., "doubling_time_hours" or "mutational_burden")
# - age_variable: the name of the age-related variable to correlate with ANOVA residuals (e.g., "age_prediction" or "age_at_sampling")
# - min_samples: the minimum number of samples required per category (default = 10)
# returns:
# - A list containing ANOVA results, residuals, correlation coefficient, p-value, and the filtered dataset.

anova_exp_settings_analysis <- function(cl_samples, response_variable, age_variable, min_samples = 10) {
  
  # check if required columns exist
  required_columns <- c(response_variable, "cancer_type", "tissue", "msi_status", age_variable)
  missing_columns <- setdiff(required_columns, colnames(cl_samples))
  if (length(missing_columns) > 0) {
    stop(paste("Error: The following columns are missing in the dataset:", paste(missing_columns, collapse = ", ")))
  }
  
  # filter dataset: keep only cancer types/tissues with at least 'min_samples' cell lines
  cl_samples <- filter_tissue_or_cancer(cl_samples, "cancer_type", "tissue", "CAccession", min_samples)
  
  # update NA values with "others"
  cl_samples$tissue_or_cancer_chosen[is.na(cl_samples$tissue_or_cancer_chosen)] <- "others"
  
  # create a status column for summary plot
  cl_samples$status <- ifelse(cl_samples$tissue_or_cancer_chosen == "others", "Others",
                              ifelse(cl_samples$tissue_or_cancer_chosen %in% cl_samples$cancer_type, "Cancer Type", "Tissue Type"))
  
  # create dynamic formula for ANOVA
  anova_formula <- as.formula(paste(response_variable, "~ as.factor(tissue_or_cancer_chosen) + as.factor(msi_status)"))
  
  # perform ANOVA
  anova_model <- aov(anova_formula, data = cl_samples)
  anova_results <- summary(anova_model)
  
  # extract residuals from ANOVA model
  cl_samples$residuals <- residuals(anova_model)
  
  # apply filtering of residuals only if response_variable is "doubling_time_hours"
  if (response_variable == "doubling_time_hours") {
    cl_samples <- cl_samples[cl_samples$residuals <= 100, ]
  }
  
  # correlation between ANOVA residuals and the selected age-related variable
  correlation <- cor(cl_samples$residuals, cl_samples[[age_variable]], use = "complete.obs")
  p_value <- cor.test(cl_samples$residuals, cl_samples[[age_variable]], method = "pearson")$p.value
  
  # print results
  print("ANOVA Results:")
  print(anova_results)
  print(paste("Correlation between ANOVA residuals and", age_variable, ":", round(correlation, 3)))
  print(paste("P-value:", ifelse(p_value < 2.2e-16, "< 2.2e-16", round(p_value, 5))))
  
  # return results as a list
  return(list(
    anova_results = anova_results,
    residuals = cl_samples$residuals,
    correlation = correlation,
    p_value = p_value,
    cl_samples_filtered = cl_samples
  ))
}






############ anova_exp_settings_cancer_specific() ###################
# The 'anova_exp_settings_cancer_specific' function performs a cancer-type specific ANOVA 
# analysis on an experimental response variable (e.g., doubling time, mutational burden) while adjusting for MSI status.
# It then calculates the correlation between ANOVA residuals and a selected age-related variable (e.g., predicted age, age at sampling).
#
# Arguments:
# - cl_samples: dataframe with required columns
# - response_variable: the experimental variable of interest (e.g., "doubling_time_hours" or "mutational_burden")
# - age_variable: the numeric variable to correlate with residuals (e.g., "age_prediction" or "age_at_sampling")
# - min_samples: minimum number of samples required per cancer type (default = 10)
#
# Returns:
# - A list containing ANOVA results, residuals, correlation coefficients, p-values, and significant findings.

anova_exp_settings_cancer_specific <- function(cl_samples, response_variable, age_variable, min_samples = 10) {
  
  set.seed(123) 
  
  # check required columns
  required_columns <- c(response_variable, "cancer_type", "msi_status", age_variable)
  missing_columns <- setdiff(required_columns, colnames(cl_samples))
  if (length(missing_columns) > 0) {
    stop(paste("Error: Missing columns in dataset:", paste(missing_columns, collapse = ", ")))
  }
  
  # filter cancer types with at least 'min_samples' samples
  cl_samples$cancer_specific <- "Excluded"
  cl_samples$cancer_specific[cl_samples$status == "Cancer Type"] <- "Included as Cancer Type"
  
  # keep only included cancer types
  cl_samples_filtered <- cl_samples[cl_samples$cancer_specific == "Included as Cancer Type", ]
  cl_samples_filtered$residuals <- NULL
  
  # unique cancer types
  unique_cancer_types <- unique(cl_samples_filtered$cancer_type)
  
  # initialize lists and results data frame
  anova_results_list <- list()
  correlation_results_list <- list()
  cumulative_correlation_results <- data.frame()
  significant_results <- data.frame(cancer_type = character(), correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)
  
  # loop through each cancer type
  for (cancer_type in unique_cancer_types) {
    
    # filter data for the current cancer type
    cancer_data <- cl_samples_filtered[cl_samples_filtered$cancer_type == cancer_type, ]
    
    # handle missing MSI status (treat NA as MSS)
    cancer_data$msi_status[is.na(cancer_data$msi_status)] <- "MSS"
    
    # ensure at least 2 samples for each MSI status level
    msi_count <- sum(cancer_data$msi_status == "MSI")
    mss_count <- sum(cancer_data$msi_status == "MSS")
    
    if (msi_count < 2) {
      cancer_data$msi_status[cancer_data$msi_status == "MSI"] <- "MSS"
    } else if (mss_count < 2) {
      cancer_data$msi_status[cancer_data$msi_status == "MSS"] <- "MSI"
    }
    
    # check if msi_status has at least two levels
    if (length(unique(cancer_data$msi_status)) > 1) {
      
      # create dynamic formula for ANOVA
      anova_formula <- as.formula(paste(response_variable, "~ as.factor(msi_status)"))
      
      # perform ANOVA
      anova_model <- aov(anova_formula, data = cancer_data)
      
      # extract residuals from ANOVA model
      cancer_data$residuals <- residuals(anova_model)
      
      # calculate correlation between residuals and selected age variable
      correlation <- cor(cancer_data$residuals, cancer_data[[age_variable]], use = "complete.obs")
      p_value <- cor.test(cancer_data$residuals, cancer_data[[age_variable]], method = "pearson")$p.value
      
      # store results
      anova_results_list[[cancer_type]] <- list(anova_model = anova_model, correlation = correlation)
      
    } else {
      # if MSI status does not have enough levels, perform direct correlation
      correlation <- cor(cancer_data[[response_variable]], cancer_data[[age_variable]], use = "complete.obs")
      p_value <- cor.test(cancer_data[[response_variable]], cancer_data[[age_variable]], method = "pearson")$p.value
      
      # store correlation results
      correlation_results_list[[cancer_type]] <- list(correlation = correlation, p_value = p_value)
    }
    
    # add results to cumulative data frame
    cumulative_correlation_results <- rbind(
      cumulative_correlation_results, 
      data.frame(
        cancer_type = cancer_type, 
        correlation = correlation, 
        p_value = p_value, 
        msi_status = unique(cancer_data$msi_status)[1]
      )
    )
    
    # store significant results (p-value < 0.05)
    if (p_value < 0.05) {
      significant_results <- rbind(significant_results, data.frame(
        cancer_type = cancer_type, correlation = correlation, p_value = p_value
      ))
    }
  }
  
  # print all correlation results
  cat("\nAll correlation results:\n")
  print(cumulative_correlation_results)
  
  # print significant results
  cat("\nSignificant results (p-value < 0.05):\n")
  if (nrow(significant_results) > 0) {
    print(significant_results)
  } else {
    print("No significant results found.")
  }
  
  # return results
  return(list(
    anova_results = anova_results_list,
    correlation_results = correlation_results_list,
    cumulative_correlation_results = cumulative_correlation_results,
    significant_results = significant_results,
    cl_samples_filtered = cl_samples_filtered
  ))
}





############ print_MSI_status() ###################
# the'print_MSI_status' function print a table that summarizes the MSI status among 
# different cancer types in a set of defined cell lines
print_MSI_status <- function(cl_samples_MSIvsMSS) {
  
  # ensure the required columns exist
  required_columns <- c("cancer_type", "msi_status")
  missing_columns <- setdiff(required_columns, colnames(cl_samples))
  
  if (length(missing_columns) > 0) {
    stop(paste("Error: Missing columns in dataset:", paste(missing_columns, collapse = ", ")))
  }
  
  # if they both exists
  msi_table <- table(cl_samples_MSIvsMSS$cancer_type, cl_samples_MSIvsMSS$msi_status)
  msi_table_df <- as.data.frame.matrix(msi_table)
  colnames(msi_table_df) <- c("MSI", "MSS")
  print(msi_table_df)
  
  return(msi_table_df)
}





############ test_msi_age_association_pan_cancer() ###################
# The 'test_msi_age_association_pan_cancer' function evaluates whether 
# cell lines with MSI (Microsatellite Instability) or MSS (Microsatellite Stability)
# are systematically older or younger based on a given age variable (predicted age or donor age)
# in a pan cancer setting
test_msi_age_association_pan_cancer <- function(cl_samples, age_variable) {
  # ensure required columns exist
  required_columns <- c(age_variable, "msi_status")
  missing_columns <- setdiff(required_columns, colnames(cl_samples))
  
  if (length(missing_columns) > 0) {
    stop(paste("Error: Missing columns in dataset:", paste(missing_columns, collapse = ", ")))
  }
  
  # check normality assumption using Shapiro-Wilk test
  shapiro_result <- shapiro.test(cl_samples[[age_variable]])
  is_normal <- shapiro_result$p.value > 0.05
  
  # perform the appropriate test
  if (is_normal) {
    test_result <- t.test(cl_samples[[age_variable]] ~ cl_samples$msi_status)
    test_type <- "T-test (parametric)"
  } else {
    test_result <- wilcox.test(cl_samples[[age_variable]] ~ cl_samples$msi_status)
    test_type <- "Wilcoxon rank-sum test (non-parametric)"
  }
  
  # extract results
  test_statistic <- round(test_result$statistic, 3)
  p_value <- ifelse(test_result$p.value < 2.2e-16, "< 2.2e-16", round(test_result$p.value, 5))
  
  # print summary
  cat("\n*** MSI vs MSS Age Distribution Analysis: ", age_variable, " ***\n", sep = "")
  cat("Normality test (Shapiro-Wilk) p-value:", round(shapiro_result$p.value, 5), "\n")
  cat("Selected test:", test_type, "\n")
  cat("Test statistic:", test_statistic, "\n")
  cat("P-value:", p_value, "\n")
  
  # return test results
  return(list(
    test_type = test_type,
    test_result = test_result,
    p_value = test_result$p.value
  ))
}





############ test_msi_age_association_cancer_specific() ###################
# The 'test_msi_age_association_cancer_specific' function evaluates whether  
# within each cancer type, cell lines with MSI (Microsatellite Instability) 
# or MSS (Microsatellite Stability) are systematically older or younger 
# based on a given age variable (predicted age or donor age).

test_msi_age_association_cancer_specific <- function(cl_samples, age_variable, min_samples = 3) {
  
  # ensure required columns exist
  required_columns <- c(age_variable, "msi_status", "cancer_type")
  missing_columns <- setdiff(required_columns, colnames(cl_samples))
  
  if (length(missing_columns) > 0) {
    stop(paste("Error: Missing columns in dataset:", paste(missing_columns, collapse = ", ")))
  }
  
  # unique cancer types
  cancer_types <- unique(cl_samples$cancer_type)
  
  # results storage
  test_results <- data.frame(cancer_type = character(), test_type = character(), p_value = numeric(), stringsAsFactors = FALSE)
  
  for (cancer in cancer_types) {
    
    # filter for each cancer
    cancer_data <- cl_samples %>% filter(cancer_type == cancer)
    
    # verify that there are enough samples
    if (nrow(cancer_data) < min_samples) {
      cat("Skipping", cancer, "- not enough samples\n")
      next
    }
    
    # normality test (Shapiro-Wilk)
    shapiro_test <- shapiro.test(cancer_data[[age_variable]])
    is_normal <- shapiro_test$p.value > 0.05
    
    # perform the appropriate test
    if (is_normal) {
      test_result <- t.test(cancer_data[[age_variable]] ~ cancer_data$msi_status)
      test_type <- "T-test (parametric)"
    } else {
      test_result <- wilcox.test(cancer_data[[age_variable]] ~ cancer_data$msi_status, exact = FALSE)
      test_type <- "Wilcoxon rank-sum test (non-parametric)"
    }
    
    # store results
    test_results <- rbind(test_results, 
                          data.frame(cancer_type = cancer, 
                                     test_type = test_type, 
                                     p_value = round(test_result$p.value, 5)))
    
    # print results
    cat("\nCancer Type:", cancer, "\n",
        "Selected test:", test_type, "\n",
        "P-value:", round(test_result$p.value, 5), "\n")
  }
  
  return(test_results)
}




############ compute_cohens_d_cancer_specific() ###################
# The 'compute_cohens_d_cancer_specific' function calculates Cohen's d effect size 
# for the difference in a given age variable (predicted age or donor age) between 
# MSI (Microsatellite Instability) and MSS (Microsatellite Stability) groups, 
# within each cancer type.
compute_cohens_d_cancer_specific <- function(cl_samples, age_variable, min_samples = 3) {
  
  # ensure required columns exist
  required_columns <- c(age_variable, "msi_status", "cancer_type")
  missing_columns <- setdiff(required_columns, colnames(cl_samples))
  
  if (length(missing_columns) > 0) {
    stop(paste("Error: Missing columns in dataset:", paste(missing_columns, collapse = ", ")))
  }
  
  # unique cancer types
  cancer_types <- unique(cl_samples$cancer_type)
  
  # results storage
  effect_size_results <- data.frame(
    cancer_type = character(), 
    cohens_d = numeric(), 
    stringsAsFactors = FALSE
  )
  
  for (cancer in cancer_types) {
    
    # filter for each cancer
    cancer_data <- cl_samples %>% filter(cancer_type == cancer)
    
    # verify that there are enough samples
    if (nrow(cancer_data) < min_samples) {
      cat("Skipping", cancer, "- not enough samples\n")
      next
    }
    
    # Extract MSI and MSS groups
    msi_group <- cancer_data %>% filter(msi_status == "MSI") %>% pull(.data[[age_variable]])
    mss_group <- cancer_data %>% filter(msi_status == "MSS") %>% pull(.data[[age_variable]])
    
    # Verify both groups have at least two samples
    if (length(msi_group) < 2 || length(mss_group) < 2) {
      cat("Skipping", cancer, "- not enough samples in both MSI and MSS groups\n")
      next
    }
    
    # Compute means and standard deviations
    mean_msi <- mean(msi_group, na.rm = TRUE)
    mean_mss <- mean(mss_group, na.rm = TRUE)
    sd_pooled <- sqrt(((length(msi_group) - 1) * var(msi_group, na.rm = TRUE) + 
                         (length(mss_group) - 1) * var(mss_group, na.rm = TRUE)) / 
                        (length(msi_group) + length(mss_group) - 2))
    
    # compute Cohen's d
    cohens_d <- (mean_msi - mean_mss) / sd_pooled
    
    # store results
    effect_size_results <- rbind(effect_size_results, 
                                 data.frame(cancer_type = cancer, 
                                            cohens_d = cohens_d))
    
    # print results
    cat("\nCancer Type:", cancer, "\n",
        "Cohen's d:", round(cohens_d, 3), "\n")
  }
  
  return(effect_size_results)
}


##########################
##### IRENE UTILS ########
##########################

library(TCGAbiolinks)


# Preprocess data
preprocessing <- function(name = "", array = "", external = TRUE,
                          dir_data = "/group/iorio/Irene/epiclock/data/",
                          dir_metadata = "/group/iorio/Irene/git_epiclock/metadata/"){
  dir_output <- ifelse(external == FALSE, paste0(dir_data, name), paste0(dir_data, "External/", name))
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("ChAMP")
  message("Loading annotation")
  if(array == "450k"){
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    xReactiveProbes <- read.csv(paste0(dir_metadata, "mask_cpgs/48639-non-specific-probes-Illumina450k.csv"))
  } 
  if(array == "EPIC"){
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    xReactiveProbes <- read.table(paste0(dir_metadata, "mask_cpgs/EPIC.hg38.mask.tsv"), header = TRUE)
    xReactiveProbes <- xReactiveProbes[xReactiveProbes$MASK_general == TRUE, 1]
  }
  
  message("Create beta matrix")
  if(external){
    library(minfi)
    print(paste0(dir_data, "External/", name, "/", name, "_RAW"))
    RGset <- read.metharray.exp(paste0(dir_data, "External/", name, "/", name, "_RAW"), verbose=T, force = TRUE)
    detP <- detectionP(RGset)
    GRset <- mapToGenome(preprocessNoob(RGset, verbose = T))
    betas <- getBeta(GRset)
    NA_position <- detP[featureNames(GRset),] > 0.05
    betas[NA_position]<- NA
  } else{
    load(paste(dir_data, name, "/data_met_raw.RData", sep=""))
    betas <- as.data.frame(assay(data_met))
  }
  n_probes <- nrow(betas)
  n_samples <- ncol(betas)
  
  message("Filter beta matrix")
  # remove cg with all NAs
  rows_to_keep <- rowSums(is.na(betas)) / ncol(betas) < 1 
  betas <- betas[rows_to_keep,]
  n_probes_keep <- nrow(betas)
  # remove samples that fail >10% of probes
  columns_to_keep <- colMeans(is.na(betas)) <= 0.1 
  betas <- betas[, columns_to_keep]
  n_samples_keep <- ncol(betas)
  # remove cg with >20% NAs
  rows_to_keep <- colMeans(is.na(betas)) <= 0.2 
  betas <- betas[rows_to_keep,]
  n_probes_keep2 <- nrow(betas)
  # remove probes: crossreactive, in sexual chromosomes and not in annepic (rs)
  probes_flt <- c(xReactiveProbes,
                  ann$Name[ann$chr %in% c("chrX","chrY")],
                  rownames(betas)[!(rownames(betas) %in% ann$Name)])
  betas_flt <- betas[!(rownames(betas) %in% probes_flt),]
  # only keep cg probes
  betas_flt_cg <- betas_flt[grep("^cg", rownames(betas_flt)),]
  n_probes_flt <- nrow(betas_flt_cg)
  # imputation
  message("Start imputation")
  library(ChAMP)
  n_nas <- sum(is.na(betas_flt_cg))
  betas_imp <- champ.impute(beta = as.matrix(betas_flt_cg), pd = NULL, method = "KNN")
  # BMIQ
  message("Start BMIQ")
  betas_bmiq <- champ.norm(betas_imp, arraytype = array, cores = 4, method = "BMIQ")
  
  # Save
  message("Saving betas preprocessed")
  save(betas_bmiq, file = paste0(dir_output, "/data_met_preprocessed.RData"))
  write(c(paste("Initial number of samples:", n_samples),
          paste("Samples with >10% failed probes:", n_samples-n_samples_keep),
          paste("Final number of samples:", n_samples_keep),
          paste("Initial number of probes:", n_probes),
          paste("Probes with all NAs:", n_probes-n_probes_keep),
          paste("Probes with >20% NAs:", n_probes_keep-n_probes_keep2),
          paste("Probes filtered:", n_probes_keep2-n_probes_flt),
          paste("Final number of probes:", n_probes_flt),
          paste("Number of NAs imputed:", n_nas)),
        file = paste0(dir_output, "/info_preprocessing.txt"))
}


# CpGs shared between Illumina methylation arrays (450k, EPIC, EPICv2)
CpGshared <- function(dir_annotations = "/group/iorio/Irene/git_epiclock/metadata/Illumina_arrays/",
                      dir_file = "/group/iorio/Irene/git_epiclock/res/resources/cpgs_shared.txt"){
  # Check if the file already exists
  if (file.exists(dir_file)) {
    return(readLines(dir_file))
  } else{
    # Load annotations
    ann_epicv2 <- data.table::fread(paste0(dir_annotations, "EPIC-8v2-0_A1.csv"), header = TRUE, fill = TRUE)[c(-1,-2,-3,-4,-5,-6,-7),]
    ann_epic <- data.table::fread(paste0(dir_annotations, "infinium-methylationepic-v-1-0-b5-manifest-file.csv"), header = TRUE, fill = TRUE)[c(-1,-2,-3,-4,-5,-6,-7),]
    ann_450k <-  data.table::fread(paste0(dir_annotations, "humanmethylation450_15017482_v1-2.csv"), header = TRUE, fill = TRUE)[c(-1,-2,-3,-4,-5,-6,-7),]
    # Remove control CpGs
    cpgs_epicv2 <- ann_epicv2$Inc.[1:which(ann_epicv2$Illumina == "[Controls]")-1]
    cpgs_epic <- ann_epic$Inc.[1:which(ann_epic$Illumina == "[Controls]")-1]
    cpgs_450k <- ann_450k$Inc.[1:which(ann_450k$Illumina == "[Controls]")-1]
    # Get CpGs shared
    cpgs_shared <- intersect(intersect(cpgs_epicv2, cpgs_epic), cpgs_450k)
    return(cpgs_shared)
  }
}


# TCGA projects
TCGA_projects <- function(){
  projects <- getGDCprojects()
  projects_tcga <- grep("^TCGA", projects$id, value = TRUE)
  return(projects_tcga)
}


# TCGA samples information
TCGA_samples <- function(dir_data = "/group/iorio/Irene/epiclock/data/",
                         project = TCGA_projects()){  
  all_samples <- NULL
  for(i in project){
    if(i != "TCGA-TGCT"){
      load(paste(dir_data, i, "/data_samples.RData", sep=""))
      sample_data_sub <- as.data.frame(sample_data)
      sample_data_sub <- sample_data_sub[, c("barcode", "patient", "age_at_index", "sample_type", "gender", "race", "ethnicity")]
      sample_data_sub$project <- i
      ifelse(is.null(all_samples), all_samples <- sample_data_sub, all_samples <- rbind(all_samples, sample_data_sub))
    }
  }
  # only samples with age annotation
  all_samples_sub <- all_samples[!is.na(all_samples$age_at_index),]
  return(all_samples_sub)
}

# TCGA beta values
TCGA_Bvalues <- function(dir_data = "/group/iorio/Irene/epiclock/data/",
                         project = TCGA_projects()){
  all_data <- NULL
  for(i in project){
    if(i != "TCGA-TGCT"){
      load(paste(dir_data, i, "/data_met_preprocessed.RData", sep=""))
      ifelse(is.null(all_data), all_data <- betas_bmiq, all_data <- merge(all_data, betas_bmiq, by = "row.names"))
      if(colnames(all_data)[1] == "Row.names"){
        rownames(all_data) <- all_data$Row.names
        all_data[,1] <- NULL
      }
    }
  }
  rm(betas_bmiq)
  all_data <- t(all_data)
  return(all_data)
}

# Predict the age
AgePred_tcga <- function(dir_predictions = "/group/iorio/Irene/git_epiclock/res/predictions/",
                    name,
                    coefs,
                    bval){
  not0_coefs <- coefs[coefs[,2] != 0, 1]
  cpgs_common <- intersect(not0_coefs, rownames(bval))
  cpgs_missing <- not0_coefs[!(not0_coefs %in% cpgs_common)][-1]
  message(paste("There are", length(cpgs_missing), "CpGs missing:",
                paste(cpgs_missing, collapse = ", ")))
  pred <- t(coefs[coefs[,1] %in% cpgs_common, 2]) %*% as.matrix(bval[cpgs_common,])
  pred <- coefs[1, 2] + pred
  write.csv(t(as.data.frame(pred)), paste0(dir_predictions, name))
  return(pred)
}





