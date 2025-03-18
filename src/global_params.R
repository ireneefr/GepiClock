####################################################
###   Global Parameters                          ###
###   Author: Digilio Alessandro                 ###
###   Date: 05/03/2025                           ###
####################################################

# in this script are collected all the global parameters

# cell lines' cancer types colors
cancer_colors <- c(
  "Acute Monocytic Leukemia" = "#ff0000",
  "Acute Myeloid Leukemia" = "#ea607f",
  "B-Cell Non-Hodgkin's Lymphoma" = "#eb8fa4",
  "B-Lymphoblastic Leukemia" = "#ebd6da",
  "Biliary Tract Carcinoma" = "#e0d14c",
  "Bladder Carcinoma" = "#b5b50e",
  "Breast Carcinoma" = "#f190f2",
  "Burkitt's Lymphoma" = "#cfb0c5",
  "Cervical Carcinoma" = "#f5c1d9",
  "Chondrosarcoma" = "#e5d69c",
  "Chronic Myelogenous Leukemia" = "#ff7fbb",
  "Colorectal Carcinoma" = "#fd8e4d",
  "Endometrial Carcinoma" = "#b59aad",
  "Esophageal Carcinoma" = "#70d3ff",
  "Esophageal Squamous Cell Carcinoma" = "#65bee6",
  "Ewing's Sarcoma" = "#75ed9b",
  "Gastric Carcinoma" = "#f5d349",
  "Glioblastoma" = "#bfb282",
  "Glioma" = "#fabf8f",
  "Head and Neck Carcinoma" = "#65cb51",
  "Hepatocellular Carcinoma" = "#60497a",
  "Hodgkin's Lymphoma" = "#ccc0da",
  "Kidney Carcinoma" = "#74f3d1",
  "Melanoma" = "#fcd4ca",
  "Mesothelioma" = "#b6cef2",
  "Neuroblastoma" = "#e2a536",
  "Non-Cancerous" = "#f2f2f2",
  "Non-Small Cell Lung Carcinoma" = "#5594f2",
  "Oral Cavity Carcinoma" = "#78f361",
  "Osteosarcoma" = "#9affba",
  "Other Blood Cancers" = "#ec4c5e",
  "Other Sarcomas" = "#da9694",
  "Other Solid Cancers" = "#d9d9d9",
  "Ovarian Carcinoma" = "#6bd6cd",
  "Pancreatic Carcinoma" = "#f84498",
  "Plasma Cell Myeloma" = "#eba6b6",
  "Prostate Carcinoma" = "#aee364",
  "Rhabdomyosarcoma" = "#f3e3a4",
  "Small Cell Lung Carcinoma" = "#176ef2",
  "Squamous Cell Lung Carcinoma" = "#91b8f2",
  "T-Cell Non-Hodgkin's Lymphoma" = "#963634",
  "T-Lymphoblastic Leukemia" = "#e2b609",
  "Thyroid Gland Carcinoma" = "#c4d79b",
  "Unknown" = "#f2f2f2"
)


tissue_colors <- c(
  "Adrenal Gland" = "#a6a6a6",
  "Biliary Tract"	= "#c6b843",
  "Bladder"	= "#d9d909",
  "Bone" =	"#69d38b",
  "Breast"	= "#cb79cc",
  "Central Nervous System" = "#988d67",
  "Cervix" = "#ff0000",
  "Endometrium" =	"#9c8595",
  "Esophagus"	= "#5ba9cd",
  "Eye" =	"#6ed0d8",
  "Haematopoietic and Lymphoid" =	"#d07d80",
  "Head and Neck" =	"#50a240",
  "Kidney" =	"#64d2b5",
  "Large Intestine" =	"#cb855d",
  "Liver" =	"#60497a",
  "Lung" =	"#7096d1",
  "Ovary" =	"#639692",
  "Pancreas" =	"#d474a1",
  "Peripheral Nervous System" =	"#b88832",
  "Placenta" =	"#92cddc",
  "Prostate" =	"#adc888",
  "Skin" =	"#deb2a6",
  "Small Intestine" =	"#fabf8e",
  "Soft Tissue" =	"#d474a1",
  "Stomach" =	"#ddbb77",
  "Testis"	= "#0071c0",
  "Thyroid"	= "#76933c",
  "Unknown" =	"#d9d9d9",
  "Uterus" =	"#c00000",
  "Vulva"	= "#eb4c5e"
)
  
# acronyms for cancer types (in the cpmments you find the counterpart in TCGA)
acronym_map <- list(
  "Colorectal Carcinoma" = "COAD",  # Colon adenocarcinoma
  "Glioblastoma" = "GBM",
  "Plasma Cell Myeloma" = "MISC",  # Not explicitly in TCGA
  "Neuroblastoma" = "MISC",  # Not explicitly in TCGA
  "Non-Small Cell Lung Carcinoma" = "LUAD",  # Lung adenocarcinoma
  "T-Lymphoblastic Leukemia" = "MISC",  # Not explicitly in TCGA
  "Hodgkin's Lymphoma" = "MISC",  # Not explicitly in TCGA
  "Head and Neck Carcinoma" = "HNSC",
  "Small Cell Lung Carcinoma" = "MISC",  # Not explicitly in TCGA
  "Pancreatic Carcinoma" = "PAAD",
  "Rhabdomyosarcoma" = "SARC",
  "Other Solid Cancers" = "MISC",
  "Other Blood Cancers" = "MISC",
  "Gastric Carcinoma" = "STAD",
  "Kidney Carcinoma" = "KIRC",  # Kidney renal clear cell carcinoma
  "Non-Cancerous" = "CNTL",  # Controls
  "Squamous Cell Lung Carcinoma" = "LUSC",
  "Esophageal Carcinoma" = "ESCA",
  "Ewing's Sarcoma" = "SARC",
  "Endometrial Carcinoma" = "UCEC",
  "Acute Myeloid Leukemia" = "LAML",
  "Glioma" = "LGG",
  "Biliary Tract Carcinoma" = "CHOL",
  "Hepatocellular Carcinoma" = "LIHC",
  "Breast Carcinoma" = "BRCA",
  "Oral Cavity Carcinoma" = "HNSC",
  "B-Cell Non-Hodgkin's Lymphoma" = "DLBC",
  "Chronic Myelogenous Leukemia" = "LCML",
  "Osteosarcoma" = "SARC",
  "Melanoma" = "SKCM",
  "Esophageal Squamous Cell Carcinoma" = "ESCA",
  "Thyroid Gland Carcinoma" = "THCA",
  "Ovarian Carcinoma" = "OV",
  "Prostate Carcinoma" = "PRAD",
  "Chondrosarcoma" = "SARC",
  "Bladder Carcinoma" = "BLCA",
  "Mesothelioma" = "MESO",
  "B-Lymphoblastic Leukemia" = "MISC",  # Not explicitly in TCGA
  "Cervical Carcinoma" = "CESC",
  "Burkitt's Lymphoma" = "DLBC",
  "Other Sarcomas" = "SARC",
  "T-Cell Non-Hodgkin's Lymphoma" = "DLBC"
)






