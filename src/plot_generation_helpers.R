####################################################
###   Plot Generation Helpers                    ###
###   Author: Digilio Alessandro                 ###
###   Date: 05/03/2025                           ###
####################################################

# in the script are written all the function to generate useful plots 

library(ggplot2)

# this plot compare the performance of the clocks in terms of correlation between the predicted age and the annotated age at sampling in cell lines
plot_clocks_performances <- function(correlations) {
  custom_colors <- c("#61D04F", "#87CEEB", "#2297E6", "#003366")
  
  correlations$Clock <- factor(
    correlations$Clock, 
    levels = c("HorvathS2013", "HannumG2013", "HorvathS2018", "LevineM2018")
  )
  
  clock_performances <- ggplot(correlations, aes(x = Clock, y = Pearson_Correlation, fill = Clock)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = custom_colors) +
    labs(x = "", y = "Pearson Correlation") +
    theme_minimal(base_size = 25) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  output_path <- paste0(figures_path, "clocks_evaluation.png")
  
  ggsave(output_path, plot = clock_performances, width = 10, height = 8, dpi = 300)
  
  return(clock_performances)
}





# UMAP visualization  of methylation data
plot_umap_methylation_data <- function(methylation_data, cl_samples, label_variable, output_filename) {
  
  # ensure the selected label_variable exists in cl_samples
  if (!(label_variable %in% colnames(cl_samples))) {
    stop("Error: the specified label variable is missing in the metadata.")
  }
  
  # handle missing or empty values
  if (label_variable == "msi_status") {
    cl_samples_filtered <- cl_samples %>% filter(msi_status %in% c("MSI", "MSS"))
  } else {
    cl_samples_filtered <- cl_samples %>% filter(!is.na(.data[[label_variable]]), .data[[label_variable]] != "")
  }
  
  # subset methylation data to match filtered cell lines
  methylation_data_filtered <- methylation_data[, colnames(methylation_data) %in% cl_samples_filtered$CAccession, drop = FALSE]
  
  # check that valid data remains
  if (ncol(methylation_data_filtered) == 0) {
    stop("Error: no matching samples found after filtering.")
  }
  
  # perform UMAP dimensionality reduction
  set.seed(123) 
  umap_result <- umap(t(methylation_data_filtered))
  
  # create UMAP dataframe
  umap_df <- as.data.frame(umap_result$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df <- cbind(umap_df, label_category = cl_samples_filtered[[label_variable]])
  
  # define plot aesthetics
  plot_title <- switch(label_variable,
                       "tissue" = "UMAP of Methylation Data - Tissue of Origin",
                       "cancer_type" = "UMAP of Methylation Data - Cancer Type",
                       "msi_status" = "UMAP of Methylation Data - MSI Status")
  
  color_label <- switch(label_variable,
                        "tissue" = "Tissue of Origin",
                        "cancer_type" = "Cancer Type",
                        "msi_status" = "MSI Status")
  
  # define color scale based on global parameters
  if (label_variable == "msi_status") {
    color_scale <- scale_color_manual(values = c("MSI" = "#2296E6", "MSS" = "#61D04F"))
  } else if (label_variable == "cancer_type") {
    color_scale <- scale_color_manual(values = cancer_colors)
  } else if (label_variable == "tissue") {
    color_scale <- scale_color_manual(values = tissue_colors)
  } else {
    stop("Error: Unsupported label variable for color mapping.")
  }
  
  # generate UMAP plot
  umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = label_category)) +
    geom_point(alpha = 0.7, size = 3) +
    labs(title = plot_title, x = "UMAP 1", y = "UMAP 2", color = color_label) +
    color_scale +
    theme_minimal(base_size = 30) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      axis.text.y = element_text(size = 18),
      axis.text.x = element_text(size = 18),
      # legend.position = "none",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 14)
    )
  
  # save the plot
  ggsave(output_filename, plot = umap_plot, width = 14, height = 8, dpi = 300)
  
  return(umap_plot)
}





# this plot shows the distribution of the predicted age (from the model's CpGs) of the cell lines
plot_distribution_PredAge <- function(cl_samples) {
distribution_PredAge <- ggplot(cl_samples, aes(x = age_prediction)) +
    geom_histogram(bins = 20, fill = "steelblue", color = "black") +
    labs(
      title = "Cell Lines' Predicted Age",
      x = "Predicted Age",
      y = "Frequency"
    ) +
    theme_minimal(base_size = 25)
  
  output_path <- paste0(figures_path, "CLs_predicted_age_distribution.png")
  
  ggsave(output_path, plot = distribution_PredAge, width = 10, height = 8, dpi = 300)
  
  return(distribution_PredAge)
}





# plot sample categories (tissue/cancer/others) after GDSC-CMP merge - ANOVA DRUG PAN
plot_sample_categories_pan_drug <- function(GDSC_age, thresh_number, output_dir) {
  GDSC_age_unique <- GDSC_age %>% distinct(CAccession, .keep_all = TRUE)
  total_samples <- nrow(GDSC_age_unique)
  
  sample_categories_pan_drug <- ggplot(GDSC_age_unique %>%
                                         mutate(cancer_type = fct_reorder(cancer_type, table(cancer_type)[cancer_type])),
                                       aes(x = cancer_type, fill = status)) +
    geom_bar(position = "dodge") +
    geom_hline(yintercept = thresh_number, linetype = "dashed", color = "red") +
    theme_minimal(base_size = 15) +
    labs(x = "", y = "Number of Cell Lines", fill = "Status") +
    scale_fill_manual(values = c("Cancer Type" = "#EEC591", "Tissue Type" = "#00CDCD", "Others" = "#DF536B")) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.x = element_text(angle = 0, size = 12, hjust = 0.7),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    coord_flip()
  
  ggsave(output_dir, plot = sample_categories_pan_drug, width = 10, height = 10, dpi = 1200)
  
  return(sample_categories_pan_drug)
}




# plot included tissue samples - ANOVA DRUG PAN
plot_included_tissue_pan_drug <- function(GDSC_age_unique) {
  filtered_data <- GDSC_age_unique %>%
    filter(status %in% c("Tissue Type", "Others")) %>%
    mutate(tissue = ifelse(status == "Others", "Others", tissue)) %>%
    mutate(tissue = fct_reorder(tissue, table(tissue)[tissue]))
  
  stacked_data <- filtered_data %>%
    group_by(tissue, cancer_type) %>%
    summarise(total_cell_lines = n(), .groups = "drop")
  
  included_tissue_pan_drug <- ggplot(stacked_data, aes(x = total_cell_lines, y = tissue, fill = cancer_type)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_d(option = "plasma", name = "Cancer Type") +
    labs(x = "Number of Cell Lines", y = NULL) +
    theme_minimal(base_size = 20) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.y = element_text(size = 17, color = "black"),
      axis.text.x = element_text(size = 15, color = "black"),
      legend.position = "none",
      legend.title = element_text(size = 17),
      legend.text = element_text(size = 17)
    )
  
  output_path <- paste0(figures_path, "included_tissue_drug_pan.png")
  ggsave(output_path, plot = included_tissue_pan_drug, width = 10, height = 8, dpi = 300)
  
  return(included_tissue_pan_drug)
}





# plot to check normality assumptions of ANOVA residuals
plot_anova_LNIC50_residuals <- function(residuals_anova_IC50, output_path) {
  # Open a PDF device
  pdf(output_path, width = 10, height = 8)  # width/height are in inches

  # Arrange plots 2x2
  par(mfrow = c(2, 2))

  # Histogram
  hist(residuals_anova_IC50,
       main = "Residuals Distribution",
       xlab = "Residuals",
       breaks = 20,
       col = "lightgray",
       border = "black")

  # Q-Q plot
  qqnorm(residuals_anova_IC50)
  qqline(residuals_anova_IC50, col = "red")

  # Boxplot per group
  boxplot(residuals_anova_IC50 ~ GDSC_age$tissue_or_cancer_chosen,
          main = "Distribution of Residuals per Group",
          xlab = "",
          ylab = "Residuals",
          col = "lightblue",
          las = 2)

  # Close device
  dev.off()
}



# volcano plot of correlation results - ANOVA DRUG PAN
plot_predage_LNIC50_correlation_pan_drug <- function(correlation_results_ANOVA, output_dir) {
  volcano_pan_drug <- ggplot(correlation_results_ANOVA, aes(x = correlation, y = -log10(adj_p_value))) +
    geom_point(aes(color = case_when(adj_p_value < 0.05 ~ correlation, TRUE ~ NA_real_)), alpha = 0.5) +
    scale_color_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0, na.value = "grey") +
    theme_minimal(base_size = 20) +
    labs(
      x = "Correlation Coefficient",
      y = "-log10(FDR-adjusted p-value)",
      title = "Volcano Plot: Correlation between Residuals and Predicted Age",
      color = "Correlation"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

  output_path <- output_dir
  ggsave(output_path, plot = volcano_pan_drug, width = 10, height = 8, dpi = 300)
  
  return(volcano_pan_drug)
}





# plot the output of the pancancer dsea in a dotplot customized to have young and old specific colors
dsea_pancancer_dotplot <- function(dsea_result, title, output_dir) {
  dsea_results_df <- as.data.frame(dsea_result)
  
  # select top 10 pathways with strongest absolute NES
  top_terms <- dsea_results_df %>%
    arrange(desc(abs(NES))) %>%
    slice(1:10) %>%
    mutate(
      Young_Old = ifelse(NES > 0, "Positive", "Negative"),
      Significance = -log10(padj)
    )
  
  # order pathways by NES (same as reorder(pathway, NES))
  ordered_targets <- top_terms %>%
    arrange(NES) %>%
    pull(pathway)
  
  color_palette <- c("Negative" = "#4C72B0", "Positive" = "#E69F00")
  
  # dotplot
  dsea_pan_dotplot <- ggplot(top_terms, aes(x = NES, y = reorder(pathway, NES), size = Significance, color = Young_Old)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = color_palette, name = "cor(pred_age, IC50)") +
    scale_size_continuous(range = c(4, 11), name = "-log10(p.adjust)") +
    labs(
      title = NULL,
      x = "Normalized Enrichment Score (NES)",
      y = "Drug Target"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.y = element_text(size = 18),
      axis.text.x = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 10)),
      size = guide_legend(override.aes = list(size = 10))
    )
  
  # save plot
  output_path <- output_dir
  ggsave(output_path, plot = dsea_pan_dotplot, width = 10, height = 8, dpi = 600)
  
  # return ordered pathways
  return(ordered_targets)
}




# stacked bar plot showing the number of drugs associated.
plot_number_target_drug <- function(significant_results_DSEA_df, dsea_pan_cancer_results, ordered_targets = NULL, output_dir, title = NULL) {
  # --- Extract and count drugs per pathway ---
  pathway_drugs <- significant_results_DSEA_df %>%
    filter(pathway %in% ordered_targets) %>%
    mutate(Drug = strsplit(leadingEdge, ", ")) %>%
    unnest(Drug) %>%
    mutate(Drug = trimws(Drug))
  
  pathway_drugs_summary <- pathway_drugs %>%
    group_by(pathway) %>%
    summarise(Drug_Count = n(), .groups = "drop")
  
  # --- Add NES info from dsea_pan_cancer_results to classify Young vs Old ---
  nes_info <- dsea_pan_cancer_results %>%
    select(pathway, NES) %>%
    distinct()
  
  pathway_drugs_summary <- pathway_drugs_summary %>%
    left_join(nes_info, by = "pathway") %>%
    mutate(
      Young_Old = ifelse(NES > 0, "Positive", "Negative"),
      pathway = factor(pathway, levels = ordered_targets)
    )
  
  # --- Colors for Young vs Old ---
  color_palette <- c("Negative" = "#4C72B0", "Positive" = "#E69F00")
  
  # --- Barplot ---
  number_drugs_plot <- ggplot(pathway_drugs_summary, aes(x = pathway, y = Drug_Count, fill = Young_Old)) +
    geom_bar(stat = "identity", width = 0.8, alpha = 0.8) +
    scale_fill_manual(values = color_palette, name = "cor(pred_age, IC50)") +
    labs(
      title = title,
      x = "Significant Drug Target",
      y = "Number of Drugs"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14)
    ) +
    coord_flip()
  
  # --- Save ---
  ggsave(output_dir, plot = number_drugs_plot, width = 12, height = 8, dpi = 600)
  
  return(number_drugs_plot)
}





# stacked bar plot showing the number and the name of drugs associated with each statistically significant drug target.
plot_target_drug_distribution <- function(significant_results_DSEA_df, ordered_targets = NULL, output_dir) {
  # extract drugs from leadingEdge
  pathway_drugs <- significant_results_DSEA_df %>%
    filter(pathway %in% ordered_targets) %>%
    mutate(Drug = strsplit(leadingEdge, ", ")) %>%
    unnest(Drug) %>%
    mutate(Drug = trimws(Drug))
  
  # count number of drugs per pathway
  pathway_drugs_summary <- pathway_drugs %>%
    group_by(pathway) %>%
    summarise(Drug_Count = n(), .groups = "drop")
  
  # join and apply pathway order
  pathway_drugs <- pathway_drugs %>%
    left_join(pathway_drugs_summary, by = "pathway") %>%
    mutate(pathway = factor(pathway, levels = ordered_targets))
  
  # stacked barplot
  drugs_involved <- ggplot(pathway_drugs, aes(x = pathway, fill = Drug)) +
    geom_bar(stat = "count", width = 0.8) +
    labs(
      title = NULL,
      x = "Significant Drug Target",
      y = "Number of Drugs",
      fill = "Drug"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      axis.text.x = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14)
    ) +
    scale_fill_viridis_d(option = "plasma") +
    coord_flip()
  
  # save plot
  output_path <- output_dir
  ggsave(output_path, plot = drugs_involved, width = 12, height = 8, dpi = 600)
  
  return(drugs_involved)
}





# bar plot showing the number of unique cancer-specific cell lines, grouped by MSI status - ANOVA DRUG CANCER SPECIFIC
plot_sample_categories_cancer_specific_drug <- function(GDSC_age_included, output_dir) {
  # unique cell lines
  GDSC_age_unique_included <- GDSC_age_included[!duplicated(GDSC_age_included$CAccession), ]
  
  total_samples <- aggregate(CAccession ~ cancer_type, data = GDSC_age_unique_included, FUN = length)
  colnames(total_samples)[2] <- "total"
  
  GDSC_age_unique_included$cancer_type <- factor(GDSC_age_unique_included$cancer_type, 
                                                 levels = names(sort(table(GDSC_age_unique_included$cancer_type), decreasing = TRUE)))
  
  sample_categories_cancer_specific_drug <- ggplot(GDSC_age_unique_included, aes(x = cancer_type, fill = msi_status)) +
    geom_bar(position = "dodge") +
    geom_text(data = total_samples, aes(x = cancer_type, y = total, label = total), 
              vjust = -0.1, size = 3, inherit.aes = FALSE) +
    theme_minimal(base_size = 15) +
    labs(x = "", y = "Number of Cell Lines", fill = "MSI Status") +
    scale_fill_manual(values = c("MSI" = "#CD8C95", "MSS" = "#6E8B3D")) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.x = element_text(angle = 0, size = 12, hjust = 0.7),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    coord_flip()
  
  ggsave(output_dir, plot = sample_categories_cancer_specific_drug, width = 10, height = 8, dpi = 1200)
  
  return(sample_categories_cancer_specific_drug)
}




# volcano plot for cumulative correlation between drug response and predicted age across cancer types
plot_predage_LNIC50_cumulative_correlation_cancer_specific_drug <- function(cumulative_correlation_results, output_dir) {
  drug_volcano_plot <- ggplot(cumulative_correlation_results, aes(x = correlation, y = -log10(adj_p_value), color = cancer_type)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 1) +
    theme_minimal(base_size = 40) +
    labs(
      x = "Correlation Coefficient",
      y = "-log10(FDR-adjusted p-value)"
    ) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.x = element_text(size = 22, face = "bold"),
      axis.text.y = element_text(size = 22, face = "bold"),
      panel.grid.major = element_line(color = "gray85"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    scale_color_manual(values = cancer_colors)
  
  ggsave(output_dir, plot = drug_volcano_plot, width = 10, height = 10, dpi = 600)
  
  return(drug_volcano_plot)
}





# volcano plot for selected cancer types 
plot_cancer_specific_volcano <- function(cumulative_correlation_results, selected_cancer_types, output_dir) {
  filtered_data <- cumulative_correlation_results %>%
    filter(cancer_type %in% selected_cancer_types)
  
  volcano_plot <- ggplot(filtered_data, aes(x = correlation, y = -log10(adj_p_value), color = cancer_type)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 1) +
    theme_minimal(base_size = 40) +
    labs(
      x = "Correlation Coefficient",
      y = "-log10(FDR-adjusted p-value)"
    ) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.x = element_text(size = 22, face = "bold"),
      axis.text.y = element_text(size = 22, face = "bold"),
      panel.grid.major = element_line(color = "gray85"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    scale_color_manual(values = cancer_colors)
  
  ggsave(
    output_dir,
    plot = volcano_plot,
    width = 10,
    height = 10,
    dpi = 600
  )
  
  return(volcano_plot)
}





# scatter plots for significant drugs found in each cancer type from the ANOVA/correlation analysis.
# each plot corresponds to a cell line for the analyzed cancer type.
plot_significant_drugs <- function(GDSC_age_filtered, significant_drugs_list, output_dir) {
  
  for (cancer_type in names(significant_drugs_list)) {
    
    # filter data for the specific cancer type
    cancer_data <- GDSC_age_filtered[GDSC_age_filtered$cancer_type == cancer_type, ]
    
    for (drug in names(significant_drugs_list[[cancer_type]])) {
      
      # filter for the specific drug
      plot_data <- subset(cancer_data, DRUG_NAME == drug)
      
      if (nrow(plot_data) > 0) {
        
        # get correlation and adj_p_value from significant_drugs_list
        correlation_value <- significant_drugs_list[[cancer_type]][[drug]]$correlation
        adj_p_value <- significant_drugs_list[[cancer_type]][[drug]]$adj_p_value
        
        # retrieve color for the specific cancer type (defined globally)
        cancer_color <- cancer_colors[[cancer_type]]
        
        # convert MIN_CONC and MAX_CONC to logarithmic scale
        min_conc_log <- log(plot_data$MIN_CONC, base = exp(1))
        max_conc_log <- log(plot_data$MAX_CONC, base = exp(1))
        
        # compute global min and max for concentrations
        min_conc_global <- min(min_conc_log, na.rm = TRUE)
        max_conc_global <- max(max_conc_log, na.rm = TRUE)
        
        # max and min LN_IC50
        min_ic50 <- min(plot_data$LN_IC50, na.rm = TRUE)
        max_ic50 <- max(plot_data$LN_IC50, na.rm = TRUE)
        
        # find the least dense region for annotation placement
        coords <- find_empty_space(plot_data$age_prediction, plot_data$LN_IC50, n_grid = 5)
        
        # create scatterplot
        drug_plot <- ggplot(plot_data, aes(x = age_prediction, y = LN_IC50)) +
          geom_point(size = 3, alpha = 0.8, color = cancer_color) +
          geom_smooth(method = "lm", color = "black") +
          geom_hline(yintercept = max_conc_global, linetype = "dashed", color = "blue", size = 1) +
          geom_hline(yintercept = min_conc_global, linetype = "dashed", color = "blue", size = 1) +
          theme_minimal(base_size = 30) +
          labs(
            title = paste(drug),
            x = "Predicted Age",
            y = "LN_IC50"
          ) +
          annotate(
            geom = "text",
            label = paste0(
              "N = ", nrow(plot_data),
              "\nR = ", round(correlation_value, 3),
              "\nadj p.val = ", ifelse(adj_p_value < 2.2e-16, "< 2.2e-16", round(adj_p_value, 3))
            ),
            x = coords["x"],
            y = coords["y"],
            size = 10,
            hjust = 0
          ) +
          theme(
            legend.position = "none",
            panel.border = element_rect(color = "black", fill = NA, size = 1.5)
          )
        
        # file path
file_path <- file.path(output_dir, paste0("Figure4C_drug_", gsub(" ", "_", drug), "_", gsub(" ", "_", cancer_type), ".pdf"))
print(file_path)  # This will print the file path, check if it's correct.

# save 
ggsave(
  filename = file_path,
  plot = drug_plot,
  width = 10,
  height = 8,
  dpi = 600
)
        
      } else {
        print(paste("No data available for", drug, "in", cancer_type))
      }
    }
  }
}





# dotplot for cancer specific DSEA results, using clustering from MutExMatSorting and coloring pathways  
# based on enrichment in young or old predicted groups.  
dsea_cancerspecific_dotplot <- function(significant_results, title, output_dir) {   
  library(MutExMatSorting)   
  library(reshape2)       
  
  color_palette <- c("Negative" = "#4C72B0", "Positive" = "#E69F00")     
  
  # classification column   
  significant_results <- significant_results %>%     
    mutate(
      NES_category = ifelse(NES > 0, "Positive", "Negative"),
      NES_category = trimws(NES_category) # remove spaces 
    )  
  
  # print unique values of NES_category to check
  print("Unique NES categories:")
  print(unique(significant_results$NES_category))
  
  # reshape the data for clustering   
  nes_matrix <- acast(significant_results, pathway ~ acronyms_cancer_type, value.var = "NES", fill = 0)   
  nes_matrix <- abs(nes_matrix)  # make all values absolute    
  
  # apply clustering using MutExMatSorting:
  # The MutExMatSorting algorithm was used to reorder pathways (rows) and cancer types (columns) 
  # in a way that minimizes overlap among signals. This highlights patterns of mutual exclusivity 
  # or co-occurrence across cancer types, making it easier to identify enriched pathways that 
  # are specific to certain cancer types. This heuristic approach optimizes the visualization 
  # by grouping similar patterns together without forcing predefined clusters. 
  sort_result <- MExMaS.MEMo(     
    nes_matrix,     
    display = TRUE,     
    cluster_cols = FALSE,     
    legend = TRUE,     
    show_rownames = TRUE,     
    show_colnames = TRUE,     
    col = c('white', 'red')   
  )       
  
  # extract ordered row and column names   
  row_order <- rev(rownames(sort_result))   
  col_order <- colnames(sort_result)       
  
  # update factor levels to apply sorting   
  significant_results$pathway <- factor(significant_results$pathway, levels = row_order)   
  significant_results$acronyms_cancer_type <- factor(significant_results$acronyms_cancer_type, levels = col_order)       
  
  # sort data accordingly   
  significant_results <- significant_results[order(significant_results$pathway, significant_results$acronyms_cancer_type), ]       
  
  # dotplot   
  dsea_dotplot <- ggplot(significant_results, aes(x = acronyms_cancer_type, y = pathway, size = abs(NES), color = NES_category)) +     
    geom_point(alpha = 0.8) +     
    scale_size(range = c(3, 10), name = "NES (absolute)") +     
    scale_color_manual(
      values = color_palette, 
      name = "cor(pred_age, IC50)", 
      breaks = c("Negative", "Positive")
    ) +     
    theme_minimal() +     
    labs(       
      title = title,       
      x = NULL,       
      y = NULL     
    ) +     
    theme_minimal(base_size = 20) +     
    theme(       
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),       
      axis.text.x = element_text(angle = 45, size = 16, hjust = 0.7),       
      axis.text.y = element_text(size = 16),       
      legend.title = element_text(size = 16),       
      legend.text = element_text(size = 16)     
    ) +     
    guides(color = guide_legend(override.aes = list(size = 8)))       
  
  ggsave(output_dir, plot = dsea_dotplot, width = 12, height = 8, dpi = 600)       
  
  return(dsea_dotplot) 
}





# bar plot showing the distribution of cell lines categorized as Cancer Type, 
# Tissue Type, or Others, with a threshold reference line.
plot_sample_categories_pan_gene <- function(depmap_combined, thresh_number, output_dir) {   
  # consider only unique cell lines
  depmap_unique <- depmap_combined %>% distinct(BROAD_ID, .keep_all = TRUE)   
  
  # generate plot
  sample_categories_pan_gene <- ggplot(depmap_unique %>% 
                                         mutate(cancer_type = fct_reorder(cancer_type, table(cancer_type)[cancer_type])),
                                       aes(x = cancer_type, fill = status)) + 
    geom_bar(position = "dodge") + 
    geom_hline(yintercept = thresh_number, linetype = "dashed", color = "red") + 
    theme_minimal(base_size = 15) + 
    labs(x = "", y = "Number of Cell Lines", fill = "Status") + 
    scale_fill_manual(values = c("Cancer Type" = "#EEC591", "Tissue Type" = "#00CDCD", "Others" = "#DF536B")) + 
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.x = element_text(angle = 0, size = 12, hjust = 0.7),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) + 
    coord_flip()
  
  ggsave(output_dir, plot = sample_categories_pan_gene, width = 10, height = 10, dpi = 1200)
  
  return(sample_categories_pan_gene)
}




# stacked bar plot  showing the number of cell lines selected for tissue analysis, 
# categorized by cancer type.
plot_included_tissue_pan_gene <- function(depmap_combined) {   
  # consider only unique cell lines
  depmap_unique <- depmap_combined %>% distinct(BROAD_ID, .keep_all = TRUE)   
  
  # filter only "Tissue Type" and "Others"
  filtered_data <- depmap_unique %>% 
    filter(status %in% c("Tissue Type", "Others")) %>% 
    mutate(tissue = ifelse(status == "Others", "Others", tissue)) %>%  # Rename 'others' as a tissue 
    mutate(tissue = fct_reorder(tissue, table(tissue)[tissue])) 
  
  # compute total cell lines per tissue-cancer combination
  stacked_data <- filtered_data %>% 
    group_by(tissue, cancer_type) %>% 
    summarise(total_cell_lines = n(), .groups = "drop") 
  
  # stacked bar plot
  included_tissue_pan_gene <- ggplot(stacked_data, aes(x = total_cell_lines, y = tissue, fill = cancer_type)) + 
    geom_bar(stat = "identity") + 
    scale_fill_viridis_d(option = "plasma", name = "Cancer Type") + 
    labs(x = "Number of Cell Lines", y = NULL) + 
    theme_minimal(base_size = 20) + 
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.y = element_text(size = 17, color = "black"),
      axis.text.x = element_text(size = 15, color = "black"),
      legend.position = "none",
      legend.title = element_text(size = 17),
      legend.text = element_text(size = 17)
    )
  
  output_path <- paste0(figures_path, "included_tissue_gene_pan.png")
  ggsave(output_path, plot = included_tissue_pan_gene, width = 10, height = 8, dpi = 300)
  
  return(included_tissue_pan_gene)
}





# volcano plot showing the correlation between ANOVA residuals of gene effect and predicted age (pan cancer)
plot_predage_gene_effect_correlation_pan_gene <- function(correlation_results_ANOVA) {   
  volcano_pan_gene <- ggplot(correlation_results_ANOVA, aes(x = correlation, y = -log10(adj_p_value))) + 
    geom_point(aes(color = case_when(adj_p_value < 0.05 ~ correlation, TRUE ~ NA_real_)), alpha = 0.5) + 
    scale_color_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0, na.value = "grey") + 
    theme_minimal(base_size = 20) + 
    labs(
      x = "Correlation Coefficient",
      y = "-log10(FDR-adjusted p-value)",
      title = "Volcano Plot: Correlation between Residuals and Predicted Age",
      color = "Correlation"
    ) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")
  
  output_path <- paste0(figures_path, "volcano_ANOVA_gene_pan_cancer.png")
  ggsave(output_path, plot = volcano_pan_gene, width = 10, height = 8, dpi = 300)
  
  return(volcano_pan_gene)
}





# plot the output of the gsea in a dotplot customized to have young and old specific colours
GSEA_yo_dotplot <- function(gsea_result, title, output_filename) {
  # convert in df
  gsea_results_df <- as.data.frame(gsea_result)
  
  # df for plot
  gsea_results_df <- gsea_results_df %>%
    mutate(
      Young_Old = ifelse(NES > 0, "Positive", "Negative"), 
      Significance = -log10(p.adjust) 
    ) %>%
    arrange(desc(abs(NES))) %>% # order per nes
    slice(1:20) %>% # Stop 20
    arrange(NES) 
  
  # how to visualize them
  gsea_results_df$Description <- factor(gsea_results_df$Description, levels = gsea_results_df$Description[order(-gsea_results_df$NES)])
  
  # colours
  color_palette <- c("Negative" = "#4C72B0", "Positive" = "#E69F00")
  
  # dotplot
  gsea_dotplot <- ggplot(gsea_results_df, aes(x = NES, y = Description, size = Significance, color = Young_Old)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = color_palette, name = "cor(pred_age, IC50") +
    scale_size_continuous(range = c(4, 11), name = "-log10(p.adjust)") +
    labs(
      title = title,
      x = "Normalized Enrichment Score (NES)",
      y = "Pathway"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.y = element_text(size = 13),
      axis.text.x = element_text(size = 13),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 10)),
      size = guide_legend(override.aes = list(size = 10))
    )  
  
  output_path <- paste0(figures_path, output_filename)
  ggsave(output_path, plot = gsea_dotplot, width = 14, height = 8, dpi = 1200)
  
  return(gsea_dotplot)
}





# plot the number of genes involved in each enrichment from GSEA in a barplot customized to have young and old specific colours
GSEA_yo_barplot <- function(gsea_result, title, output_filename) {
  # convert resulta in a df
  gsea_results_df <- as.data.frame(gsea_result)
  
  # df for the plot
  gsea_results_df <- gsea_results_df %>%
    mutate(
      Young_Old = ifelse(NES > 0, "Positive", "Negative"), 
      Significance = -log10(p.adjust)
    ) %>%
    arrange(desc(abs(NES))) %>% # order per NES
    slice(1:20) %>% # top 20
    arrange(NES) 
  
  # # how to visualize them
  gsea_results_df$Description <- factor(gsea_results_df$Description, levels = gsea_results_df$Description[order(-gsea_results_df$NES)])
  
  # colours
  color_palette <- c("Negative" = "#4C72B0", "Positive" = "#E69F00")
  
  # barplot
  gsea_barplot <- ggplot(gsea_results_df, aes(x = setSize, y = Description, fill = Young_Old)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = color_palette, name = "cor(pred_age, IC50)") +
    labs(
      title = title,
      x = "Number of Genes",
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.y = element_text(size = 13),
      axis.text.x = element_text(size = 13),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13)
    )
  
  output_path <- paste0(figures_path, output_filename)
  ggsave(output_path, plot = gsea_barplot, width = 14, height = 8, dpi = 1200)
  
  return(gsea_barplot)
}





# bar plot showing the number of unique cancer-specific cell lines, grouped by MSI status - ANOVA GENE CANCER SPECIFIC
plot_sample_categories_cancer_specific_gene <- function(gene_dep_age_filtered, output_dir) {
  # unique cell lines
  gene_dep_age_unique_filtered<- gene_dep_age_filtered[!duplicated(gene_dep_age_filtered$CAccession), ]
  
  total_samples <- aggregate(CAccession ~ cancer_type, data = gene_dep_age_unique_filtered, FUN = length)
  colnames(total_samples)[2] <- "total"
  
  gene_dep_age_filtered <- gene_dep_age_filtered %>%
    mutate(cancer_type = fct_reorder(cancer_type, table(cancer_type)[cancer_type], .desc = TRUE))
  
  cancer_plot <- ggplot(gene_dep_age_filtered, aes(x = cancer_type, fill = msi_status)) +
    geom_bar(position = "dodge") +
    geom_text(data = total_samples, aes(x = cancer_type, y = total, label = total),
              vjust = -0.1, size = 3, inherit.aes = FALSE) +
    theme_minimal(base_size = 15) +
    labs(x = "", y = "Number of Cell Lines", fill = "MSI Status") +
    scale_fill_manual(values = c("MSI" = "#CD8C95", "MSS" = "#6E8B3D")) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.x = element_text(angle = 0, size = 12, hjust = 0.7),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    coord_flip()
  
  ggsave(output_dir, plot = cancer_plot, width = 10, height = 8, dpi = 1200)
  
  return(cancer_plot)
}



# volcano plot for cumulative correlation between gene dependencies and predicted age across cancer types
plot_predage_geneDep_cumulative_correlation_cancer_specific <- function(cumulative_correlation_results) {  
  gene_volcano_plot <- ggplot(cumulative_correlation_results, aes(x = correlation, y = -log10(adj_p_value), color = cancer_type)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 1) +
    theme_minimal(base_size = 40) +
    labs(
      x = "Correlation Coefficient",
      y = "-log10(FDR-adjusted p-value)"
    ) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.x = element_text(size = 22, face = "bold"),
      axis.text.y = element_text(size = 22, face = "bold"),
      panel.grid.major = element_line(color = "gray85"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    scale_color_manual(values = cancer_colors)  
  
  # save the plot
  ggsave(paste0(figures_path, "volcano_cumulative_geneDep_cancer_specific.png"), 
         plot = gene_volcano_plot, width = 10, height = 10, dpi = 300)
  
  return(gene_volcano_plot)  
}





# volcano plot for selected cancer types
plot_cancer_specific_geneDep_volcano <- function(cumulative_correlation_results, selected_cancer_types, output_dir) {  
  filtered_data <- cumulative_correlation_results %>%
    filter(cancer_type %in% selected_cancer_types)  
  
  volcano_plot <- ggplot(filtered_data, aes(x = correlation, y = -log10(adj_p_value), color = cancer_type)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 1) +
    theme_minimal(base_size = 40) +
    labs(
      x = "Correlation Coefficient",
      y = "-log10(FDR-adjusted p-value)"
    ) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.x = element_text(size = 22, face = "bold"),
      axis.text.y = element_text(size = 22, face = "bold"),
      panel.grid.major = element_line(color = "gray85"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    scale_color_manual(values = cancer_colors)  
  
  # save the plot
  ggsave(output_dir, plot = volcano_plot, width = 10, height = 10, dpi = 600)
  
  return(volcano_plot)  
}





# scatter plots for significant genes found in each cancer type from the anova/correlation analysis.
# each plot corresponds to a cell line for the analyzed cancer type.
plot_significant_genes <- function(gene_dep_age_filtered, significant_genes_list, output_dir) {
  
  for (cancer_type in names(significant_genes_list)) {
    
    # filter data for the specific cancer type
    cancer_data <- gene_dep_age_filtered[gene_dep_age_filtered$cancer_type == cancer_type, ]
    
    for (gene in names(significant_genes_list[[cancer_type]])) {
      
      # filter for the specific gene
      plot_data <- subset(cancer_data, Gene == gene)
      
      if (nrow(plot_data) > 0) {
        
        # get correlation and adj_p_value from significant_genes_list
        correlation_value <- significant_genes_list[[cancer_type]][[gene]]$correlation
        adj_p_value <- significant_genes_list[[cancer_type]][[gene]]$adj_p_value
        
        # retrieve color for the specific cancer type
        cancer_color <- cancer_colors[[cancer_type]]
        
        # find the least dense region for annotation placement
        coords <- find_empty_space(plot_data$age_prediction, plot_data$gene_effect, n_grid = 5)
        
        # scatter plot
        gene_plot <- ggplot(plot_data, aes(x = age_prediction, y = gene_effect)) +
          geom_point(size = 3, alpha = 0.8, color = cancer_color) +
          geom_smooth(method = "lm", color = "black") +
          theme_minimal(base_size = 30) +
          labs(
            title = paste(gene),
            x = "Predicted Age",
            y = "Gene Dependency Score"
          ) +
          annotate(
            geom = "text",
            label = paste0(
              "N = ", nrow(plot_data),
              "\nR = ", round(correlation_value, 3),
              "\nadj p.val = ", ifelse(adj_p_value < 2.2e-16, "< 2.2e-16", round(adj_p_value, 5))
            ),
            x = coords["x"],
            y = coords["y"],
            size = 10,
            hjust = 0
          ) +
          theme(
            legend.position = "none",
            panel.border = element_rect(color = "black", fill = NA, size = 1.5)
          )
        
        # save 
        ggsave(
          paste0(output_dir, "Figure5B_gene_", gsub(" ", "_", gene), "_", gsub(" ", "_", cancer_type), ".pdf"),
          plot = gene_plot,
          width = 10,
          height = 8,
          dpi = 600
        )
        
      } else {
        print(paste("No data available for", gene, "in", cancer_type))
      }
    }
  }
}





# dotplot for cancer specific GSEA (top 3 pathways per cancer type), using clustering from MutExMatSorting and coloring pathways  
# based on enrichment in young or old predicted groups.  
gsea_cancerspecific_dotplot <- function(significant_results_bp_df, title, output_dir) {
  
  library(MutExMatSorting)
  library(reshape2)
  
  set.seed(123) 
  
  color_palette <- c("Negative" = "#4C72B0", "Positive" = "#E69F00")
  
  # classify pathways based on NES sign
  significant_results <- significant_results_bp_df %>%
    filter(p.adjust < 0.05) %>%  # keep only significant pathways
    mutate(
      NES_category = ifelse(NES > 0, "Positive", "Negative"),
      NES_category = trimws(NES_category) 
    ) %>%
    group_by(cancer_type) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = 3) %>%  # select top 3 pathways per cancer type
    ungroup()

  # wrap pathway names to max 4 words
  significant_results <- significant_results %>%
  mutate(Description = sapply(strsplit(Description, " "), function(words) {
    if (length(words) > 4) {
        # insert a line break every 4 words
        paste(sapply(seq(1, length(words), by = 4), function(i)
        paste(words[i:min(i + 3, length(words))], collapse = " ")),
            collapse = "\n"
        )
        } else {
        paste(words, collapse = " ")
        }
    }))

  # print unique NES categories to check
  print("Unique NES categories:")
  print(unique(significant_results$NES_category))
  
  # reshape data for clustering
  nes_matrix <- acast(significant_results, Description ~ cancer_type, value.var = "NES", fill = 0)
  nes_matrix <- abs(nes_matrix)  # take absolute NES for clustering
  
  # apply MutExMatSorting clustering
  sort_result <- MExMaS.MEMo(
    nes_matrix,
    display = TRUE,
    cluster_cols = FALSE,
    legend = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    col = c('white', 'red')
  )
  
  # extract ordered row (pathway) and column (cancer type) names
  row_order <- rev(rownames(sort_result))
  col_order <- colnames(sort_result)
  
  # apply sorting to the data
  significant_results$Description <- factor(significant_results$Description, levels = row_order)
  significant_results$cancer_type <- factor(significant_results$cancer_type, levels = col_order)
  
  # sort data accordingly
  significant_results <- significant_results[order(significant_results$Description, significant_results$cancer_type), ]
  
  # dotplot
  gsea_dotplot <- ggplot(significant_results, aes(x = cancer_type, y = Description, size = setSize, color = NES_category)) +
    geom_point(alpha = 0.8) +
    scale_size(range = c(3, 10), name = "Gene Ratio") +
    scale_color_manual(
      values = color_palette,
      name = "cor(pred_age, IC50)",
      breaks = c("Negative", "Positive")
    ) +
    theme_minimal() +
    labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 20) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.x = element_text(angle = 45, size = 11, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 11)
    ) +
    guides(color = guide_legend(override.aes = list(size = 8)))
  
  # save the plot
  ggsave(output_dir, plot = gsea_dotplot, width = 14, height = 13, dpi = 600)
  
  return(gsea_dotplot)
}















# plot that shows landscape comparison among cancer types, showing absolute Pearson correlation values 
# for a given experimental setting variable (e.g., doubling time, mutational burden) against predicted age and donor age.
plot_landscape_correlation_exp_settings_age <- function(cl_samples_filtered, response_variable, age_variable_pred, age_variable_sample, cancer_colors) {
  
  # compute Absolute Correlations for Each Cancer Type
  correlation_data <- data.frame(
    cancer_type = unique(cl_samples_filtered$cancer_type),
    abs_corr_pred_age = sapply(unique(cl_samples_filtered$cancer_type), function(cancer) {
      subset_data <- cl_samples_filtered[cl_samples_filtered$cancer_type == cancer, ]
      abs(cor(subset_data[[response_variable]], subset_data[[age_variable_pred]], use = "complete.obs"))
    }),
    abs_corr_sample_age = sapply(unique(cl_samples_filtered$cancer_type), function(cancer) {
      subset_data <- cl_samples_filtered[cl_samples_filtered$cancer_type == cancer, ]
      abs(cor(subset_data[[response_variable]], subset_data[[age_variable_sample]], use = "complete.obs"))
    })
  )
  
  # count how many cancer types have a higher absolute correlation in predicted age
  higher_in_predicted_age <- sum(correlation_data$abs_corr_pred_age > correlation_data$abs_corr_sample_age)
  total_cancer_types <- nrow(correlation_data)
  
  cat("Number of cancer types where absolute correlation is higher in predicted age:", higher_in_predicted_age, "\n")
  cat("Total number of cancer types analyzed:", total_cancer_types, "\n")
  
  # define axis labels dynamically based on the response variable
  x_label <- paste("(Abs) Correlation:", response_variable, "vs", age_variable_pred)
  y_label <- paste("(Abs) Correlation:", response_variable, "vs", age_variable_sample)
  plot_filename <- paste0(figures_path, "landscape_correlation_", response_variable, ".png")
  
  # scatterplot
  correlation_plot <- ggplot(correlation_data, aes(x = abs_corr_pred_age, y = abs_corr_sample_age, color = cancer_type)) +
    geom_point(size = 4, alpha = 0.8) +  
    scale_color_manual(values = cancer_colors) +  
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +  
    scale_x_continuous(limits = c(0, 1)) +  
    scale_y_continuous(limits = c(0, 1)) +  
    labs(
      x = x_label,
      y = y_label,
      color = "Cancer Type"
    ) +
    theme_minimal(base_size = 25) +  
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.position = "none" # "right",
      # legend.title = element_text(size = 20),
      # legend.text = element_text(size = 15)
    )
  
  ggsave(filename = plot_filename, plot = correlation_plot, width = 10, height = 10, dpi = 300)
  
  return(correlation_plot)
}





# this boxplot compares the distribution of a selected age-related variable between MSI and MSS groups
boxplot_msi_age_pancancer <- function(cl_samples, age_variable, p_value = NULL) {
  
  # ensure required columns exist
  if (!(age_variable %in% colnames(cl_samples)) || !"msi_status" %in% colnames(cl_samples)) {
    stop("Error: The specified age variable or 'msi_status' column is missing in the dataset.")
  }
  
  # compute p-value if not provided
  if (is.null(p_value)) {
    test_result <- test_msi_age_association_pan_cancer(cl_samples, age_variable)
    p_value <- test_result$p_value
  }
  
  # format the p-value for display
  if (!is.na(p_value) && !is.null(p_value)) {
    p_value_text <- paste0("p-value = ", signif(p_value, 3))
  } else {
    p_value_text <- "p-value not available"
  }
  
  # boxplot
  boxplot_age_msi <- ggplot(cl_samples, aes(x = msi_status, y = .data[[age_variable]], fill = msi_status)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    theme_minimal(base_size = 25) +
    scale_fill_manual(values = c("MSI" = "#2296E6", "MSS" = "#61D04F")) +
    labs(
      title = NULL,
      subtitle = p_value_text,
      x = "MSI Status",
      y = age_variable, 
      fill = "MSI Status"
    ) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
  
  # save
  file_name <- paste0(figures_path, "MSI_pan_", age_variable, ".png")
  ggsave(filename = file_name, 
         plot = boxplot_age_msi, width = 10, height = 8, dpi = 300)
  
  return(boxplot_age_msi)
}





# plot a boxplots for each cancer type  showing the distribution of a selected age variable 
# between MSI and MSS groups
boxplot_msi_age_cancer_specific <- function(cl_samples, age_variable, output_dir = "results/plots/MSI_cancer_type/") {
  
  # ensure required columns exist
  if (!(age_variable %in% colnames(cl_samples)) || !"msi_status" %in% colnames(cl_samples) || !"cancer_type" %in% colnames(cl_samples)) {
    stop("Error: The specified age variable or 'msi_status'/'cancer_type' column is missing in the dataset.")
  }
  
  # create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # unique cancer types
  cancer_types <- unique(cl_samples$cancer_type)
  
  for (cancer in cancer_types) {
    
    # filter for each cancer
    cancer_data <- cl_samples %>% filter(cancer_type == cancer)
    
    # verify that there are enough samples
    if (nrow(cancer_data) < 3) {
      cat("Skipping", cancer, "- not enough samples\n")
      next
    }
    
    # run test to get p-value
    test_results <- test_msi_age_association_cancer_specific(cancer_data, age_variable)
    p_value <- test_results$p_value[test_results$cancer_type == cancer]
    
    # format the p-value for display
    p_value_text <- ifelse(!is.na(p_value) && !is.null(p_value), 
                           paste0("p-value = ", formatC(p_value, format = "f", digits = 5)), 
                           "p-value not available")
    
    # count number of samples
    sample_counts <- cancer_data %>%
      group_by(msi_status) %>%
      summarise(count = n(), max_y = max(.data[[age_variable]], na.rm = TRUE), .groups = "drop")
    
    y_max <- max(cancer_data[[age_variable]], na.rm = TRUE) + 3
    
    # boxplot
    p <- ggplot(cancer_data, aes(x = msi_status, y = .data[[age_variable]], fill = msi_status)) +
      geom_boxplot(alpha = 0.9, outlier.shape = NA) +  
      theme_minimal(base_size = 30) +
      scale_fill_manual(values = c("MSI" = "#2296E6", "MSS" = "#61D04F")) +  
      labs(title = cancer,
           subtitle = p_value_text,
           x = "MSI Status",
           y = age_variable) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        legend.position = "none"
      ) +
      geom_text(data = sample_counts, aes(x = msi_status, y = y_max, label = paste0("n=", count)), size = 7, vjust = 0)
    
    # save plot
    output_path <- paste0(output_dir, gsub(" ", "_", cancer), "_MSI_", age_variable, ".png")
    ggsave(output_path, plot = p, width = 8, height = 8, dpi = 300)
    
    cat("Saved plot for", cancer, "\n")
  }
}





# this function generates a lollipop plot to visualize Cohen's d for MSI vs MSS age differences
lollipop_cohens_d_cancer_specific <- function(test_results_pred_age, test_results_donor_age) {
  
  # ensure necessary columns exist in the input data
  required_columns <- c("cancer_type", "cohens_d")
  if (!all(required_columns %in% colnames(test_results_pred_age)) || 
      !all(required_columns %in% colnames(test_results_donor_age))) {
    stop("Error: The required columns 'cancer_type' and 'cohens_d' are missing in the input datasets.")
  }
  
  # combine predicted age and donor age Cohen's d into a single dataframe
  cohens_d_data <- data.frame(
    cancer_type = rep(test_results_pred_age$cancer_type, 2),
    cohens_d = c(test_results_pred_age$cohens_d, test_results_donor_age$cohens_d),
    age_type = rep(c("Predicted Age", "Age at Sampling"), each = nrow(test_results_pred_age))
  )
  
  # define colors for predicted age and donor age
  age_colors <- c("Predicted Age" = "#E89E84", "Age at Sampling" = "#E3C16F")
  
  # generate lollipop plot
  lollipop_plot <- ggplot(cohens_d_data, aes(x = reorder(cancer_type, cohens_d), y = cohens_d, color = age_type)) +
    geom_segment(aes(xend = cancer_type, y = 0, yend = cohens_d), linewidth = 1) +
    geom_point(size = 4) +
    scale_color_manual(values = age_colors) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    coord_flip() +  # flip x and y axes for readability
    theme_minimal(base_size = 20) +
    labs(title = "Cohen's d for MSI vs MSS Across Cancer Types",
         x = "Cancer Type", 
         y = "Cohen's d",
         color = "Age Variable") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))
  
  # define output file path using global figures_path
  output_file <- paste0(figures_path, "CohenD_MSI_cancer_type.png")
  
  # save plot
  ggsave(output_file, plot = lollipop_plot, width = 12, height = 8, dpi = 300)
  
  cat("Lollipop plot saved at:", output_file, "\n")
  
  return(lollipop_plot)
}








