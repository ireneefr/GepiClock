##########################################################
###   Characterization of CpGs included in the model   ###
###   Author: Irene Fernández Rebollo                  ###
###   Date: 20/11/2024                                 ###
##########################################################

setwd("/Volumes/iorio/Irene/epiclock_dev")
source("src/utils.R")
dir_metadata <- "metadata/"

# Load packages
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(ggplot2)
library(entropy)

# Load data
model_coefs <- read.csv("results/TCGA/model/model_coefs.csv")
model_coefs <- subset(model_coefs, s0 != 0) #4863  2
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_model_cpgs <- as.data.frame(merge(ann, model_coefs, by.x = "row.names", by.y = "X"))


chr_ann <- as.data.frame(table(ann[!(ann$chr %in% c("chrX", "chrY")), "chr"]))
chr_ann$Perc_ann <- chr_ann$Freq/nrow(ann[!(ann$chr %in% c("chrX", "chrY")),])*100
chr_ann_model <- as.data.frame(table(ann_model_cpgs$chr))
chr_ann_model$Perc_ann <- chr_ann_model$Freq/nrow(ann_model_cpgs)*100
df <- merge(chr_ann_model, chr_ann, by = "Var1")

# 1. Define totals
total_clock <- sum(chr_ann_model$Freq)
total_array <- sum(chr_ann$Freq)

# 2. Function to run Fisher for each row
results <- apply(df, 1, function(row) {
  # Extract numbers for this specific chromosome
  n_clock_chr <- as.numeric(row['Freq.x'])
  n_array_chr <- as.numeric(row['Freq.y'])
  
  # Construct the 2x2 table
  # [Clock_on_Chr, Clock_off_Chr]
  # [Array_not_Clock_on_Chr, Array_not_Clock_off_Chr]
  mat <- matrix(c(
    n_clock_chr,                               # In Clock, On Chr
    total_clock - n_clock_chr,                 # In Clock, Off Chr
    n_array_chr - n_clock_chr,                 # Not in Clock, On Chr
    (total_array - n_array_chr) - (total_clock - n_clock_chr) # Not in Clock, Off Chr
  ), nrow = 2)
  
  f_test <- fisher.test(mat)
  return(c(p_value = f_test$p.value, odds_ratio = f_test$estimate[[1]]))
})

# 3. Add results back to your dataframe
df$p_value <- results[1,]
df$odds_ratio <- results[2,]

# 4. Correct for multiple testing (CRITICAL in genomics)
df$p_adj <- p.adjust(df$p_value, method = "BH")
df$significance <- ifelse(df$p_adj <= 0.05, "p_adj <= 0.05", "p_adj > 0.05")

# Plot 
ggplot(df, aes(x = Perc_ann.x, y = Perc_ann.y, label = Var1, fill = significance)) +
  geom_abline() +
  geom_label() +
  scale_fill_manual(values = c('red3', 'grey80'), name = "Fisher test") +
  xlab("% CpGs/Chr in model") + ylab("% CpGs/Chr in Illumina array") +
  theme_minimal(base_size = 15)
ggsave("Supplementary Figures/Supplementary_Figure_18.pdf", width = 9, height = 7)


# Check model CpGs distribution across chromosomes
mean(df$Freq.x/sum(df$Freq.x)*100)
sd(df$Freq.x/sum(df$Freq.x)*100)
hist(df$Freq.x/sum(df$Freq.x)*100)

H <- entropy(df$Freq.x/sum(df$Freq.x), unit = "log2") #Shannon entropy (uncertainty of a distribution)
H_max <- log2(length(df$Freq.x)) #Maximum possible entropy with n categories
H_norm <- H / H_max #Normalized entropy
H_norm #close to 1 - distribution very uniform, close to 0 - distribution very concentrated

# head(ann_model_cpgs)
# chr_weights <- ann_model_cpgs %>%
#   group_by(chr) %>%
#   summarise(sum = sum(abs(s0)),
#             count = n())
# chr_weights$density <- chr_weights$sum/chr_weights$count
# ggplot(chr_weights, aes(x = chr, y = density)) +
#   geom_col()
