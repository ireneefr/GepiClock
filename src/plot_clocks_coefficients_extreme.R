setwd("/Volumes/iorio/Irene/epiclock_dev")
source("src/utils.R")
dir_results <- "results/TCGA/"
dir_metadata <- "metadata/"

res_test <- as.data.frame(read_csv(paste0(dir_results, "model/test_clock_cpgs_extreme.csv")))
res_test$data <- "Primary Tumour (Test)"
res_tumor <- as.data.frame(read_csv(paste0(dir_results, "model/tumor_clock_cpgs_extreme.csv")))
res_tumor$data <- "Primary Tumour (Others)"
res_normal <- as.data.frame(read_csv(paste0(dir_results, "model/normal_clock_cpgs_extreme.csv")))
res_normal$data <- "Normal Tissue"

# Summarize information
res <- rbind(res_test, res_normal, res_tumor)
clocks_performance <- as.data.frame(read_csv(paste0(dir_results, "model/clocks_performance.csv")))
clocks_performance <- clocks_performance[clocks_performance$Clock != "RebolloI2025",]

# Plot
color_type <- c("Primary Tumour (Test)" = "steelblue1",
                "Primary Tumour (Others)" = "steelblue3",
                "Normal Tissue" = "chartreuse3")
# png(filename = "/group/iorio/Irene/epiclock/plots/missing_coefficients.png",
#     width = 12, height = 7, units = 'in', res = 600)
ggplot(res, aes(x = Clock_name, y = mean, fill = data, color = data)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.02, position = position_dodge(width = 0.8), color = "black") +
  scale_fill_manual(values = color_type, name = "Data subset") +
  scale_color_manual(values = color_type, name = "Data subset") +
  theme_bw(base_size = 20) +
  xlab("Missing Coefficients (%)") + ylab("Pearson Correlation (R)") 
# dev.off()

png(filename = "/Volumes/iorio/Irene/epiclock/plots/rebollo_extremevsclocks_performance.png",
    width = 12, height = 7, units = 'in', res = 600)
ggplot(clocks_performance, aes(x = Clock, y = R, fill = Sample_type, color = Sample_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  scale_fill_manual(values = color_type, name = "Data subset") +
  scale_color_manual(values = color_type, name = "Data subset") +
  
  geom_point(data = res, aes(x = Clock_name, y = R, group = data, shape = "Pearson with Rebollo coefficients"), size = 3,
             position = position_dodge(width = 0.8), color = "red", inherit.aes = FALSE) +
  scale_shape_manual(name = "Extreme CpGs subset", values = 16) +
  
  scale_x_discrete(limits=rev) +
  theme_bw(base_size = 20) +
  xlab("") + ylab("Pearson Correlation (R)") +
  coord_flip()
dev.off()
