###################################
###   Supplementary Figure 15   ###
###################################

setwd("/Volumes/iorio/Irene/epiclock_dev")
source("src/utils.R")
dir_results <- "results/TCGA/"
dir_metadata <- "metadata/"

res_test <- as.data.frame(read_csv(paste0(dir_results, "model/test_clock_cpgs.csv")))
res_test$data <- "Test"
res_tumor <- as.data.frame(read_csv(paste0(dir_results, "model/tumor_clock_cpgs.csv")))
res_tumor$data <- "Tumor"
res_normal <- as.data.frame(read_csv(paste0(dir_results, "model/normal_clock_cpgs.csv")))
res_normal$data <- "Normal"

# Summarize information
res_test_summary <- as.data.frame(res_test) %>%
  group_by(Clock_name) %>%
  summarise(
    sd = as.numeric(sd(R, na.rm = TRUE)),
    mean = as.numeric(mean(R)),
    mean_mae = as.numeric(mean(MAE)),
    data = "Primary Tumour (Test)"
  )
res_normal_summary <- as.data.frame(res_normal) %>%
  group_by(Clock_name) %>%
  summarise(
    sd = as.numeric(sd(R, na.rm = TRUE)),
    mean = as.numeric(mean(R)),
    mean_mae = as.numeric(mean(MAE)),
    data = "Normal Tissue"
  )
res_tumor_summary <- as.data.frame(res_tumor) %>%
  group_by(Clock_name) %>%
  summarise(
    sd = as.numeric(sd(R, na.rm = TRUE)),
    mean = as.numeric(mean(R)),
    mean_mae = as.numeric(mean(MAE)),
    data = "Primary Tumour (Others)"
  )
res <- rbind(res_test_summary, res_normal_summary, res_tumor_summary)
clocks_performance <- as.data.frame(read_csv(paste0(dir_results, "model/clocks_performance.csv")))
clocks_performance <- clocks_performance[clocks_performance$Clock != "RebolloI2025",]

# Plot
color_type <- c("Primary Tumour (Test)" = "steelblue1",
                "Primary Tumour (Others)" = "blueviolet",
                "Normal Tissue" = "chartreuse3")

pdf("/Volumes/iorio/Irene/epiclock_dev/Supplementary Figures/Supplementary_Figure_15.pdf",
    width = 12, height = 7)
ggplot(clocks_performance, aes(x = Clock, y = R, fill = Sample_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  scale_fill_manual(values = color_type, name = "Data subset",
                    labels = c("Normal samples\n(Independent set)", "Tumour samples\n(Independent set)", "Tumour samples\n(Test set)")) +
  guides(fill = guide_legend(reverse = TRUE)) +

  geom_errorbar(data = res, aes(x = Clock_name, y = mean, ymin = mean - sd, ymax = mean + sd, group = data),
                width = 0.1, position = position_dodge(width = 0.8), color = "black", inherit.aes = FALSE) +
  geom_point(data = res, aes(x = Clock_name, y = mean, group = data, shape = "Mean ± SD"), size = 3,
             position = position_dodge(width = 0.8), color = "red", inherit.aes = FALSE) +
  scale_shape_manual(name = "Random CpGs subset", values = 16) +
  
  scale_x_discrete(limits=rev) +
  theme_bw(base_size = 20) +
  theme(axis.text = element_text(color = "black")) +
  xlab("") + ylab("Pearson Correlation (R)") +
  coord_flip()
dev.off()
