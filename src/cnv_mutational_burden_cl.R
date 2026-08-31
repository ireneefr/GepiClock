########################################################
###   Mutational/CNV burden vs Age in Cell Lines     ###
###   Author: Irene Fernández Rebollo                ###
###   Date: 25/08/2026                               ###
########################################################

setwd("/group/iorio/Irene/epiclock_dev")
source("src/utils.R")
library(caret)
library(dplyr)
library(ggpubr)

# Load data
mut_burden <- read.csv("metadata/CMP_annotations/model_annotations_predicted_age.csv")
cnv <- read.csv("metadata/cnv/OmicsGlobalSignatures.csv")

# Plot mutational load vs age
pdf("results/mutational_burden/cl_mut_burden_age.pdf", height = 8, width = 8)
ggplot(mut_burden, aes(x = age_at_sampling, y = mutational_burden)) +
  geom_point() +
  stat_cor(aes(group = 1), method = "spearman", colour = "black",
           label.x.npc = "left", label.y.npc = "top", cor.coef.name = "rho") +
  xlab("Age at sampling") +
  ylab("Mutational burden") +
  theme_minimal(base_size = 20)
dev.off()
pdf("results/mutational_burden/cl_mut_burden_gepiclock.pdf", height = 8, width = 8)
ggplot(mut_burden, aes(x = age_prediction, y = mutational_burden)) +
  geom_point() +
  stat_cor(aes(group = 1), method = "spearman", colour = "black",
           label.x.npc = "left", label.y.npc = "top", cor.coef.name = "rho") +
  xlab("GepiClock Age") +
  ylab("Mutational burden") +
  theme_minimal(base_size = 20)
dev.off()

# Plot CNV burden vs age
cnv_mut_burden <- merge(cnv, mut_burden, by.x = "ModelID", by.y = "BROAD_ID")
pdf("results/mutational_burden/cl_cnv_burden_age.pdf", height = 8, width = 8)
ggplot(cnv_mut_burden, aes(x = age_at_sampling, y = CIN)) +
  geom_point() +
  stat_cor(aes(group = 1), method = "spearman", colour = "black",
           label.x.npc = "left", label.y.npc = "top", cor.coef.name = "rho") +
  xlab("Age at sampling") +
  ylab("Mutational burden") +
  theme_minimal(base_size = 20)
dev.off()
pdf("results/mutational_burden/cl_cnv_burden_gepiclock.pdf", height = 8, width = 8)
ggplot(cnv_mut_burden, aes(x = age_prediction, y = CIN)) +
  geom_point() +
  stat_cor(aes(group = 1), method = "spearman", colour = "black",
           label.x.npc = "left", label.y.npc = "top", cor.coef.name = "rho") +
  xlab("GepiClock Age") +
  ylab("Mutational burden") +
  theme_minimal(base_size = 20)
dev.off()
