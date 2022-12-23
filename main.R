####
# Main script for Project CPIDM
####


# Configuration -----------------------------------------------------------
ver <- "1222"
setwd("D:/_Projects/CPIDM/R_analysis/")
# Creating output parent directory
if (!dir.exists(paste0("./output/", ver))) {
  dir.create(paste0("./output/", ver))
  print(paste("Output directory for version", ver, "created"))
} else {
  print("DIRECTORY ALREADY EXISTS!!!")
}


# Import Required Libraries -----------------------------------------------
source("./R_script/pre/import_total.R")


# Common Process ----------------------------------------------------------
# Reads annotation.csv and creates sampleTable object
source("./R_script/pre/common_pre.R")


# DESeq2 ------------------------------------------------------------------
# Create output directory
if (!dir.exists(paste0("./output/", ver, "/MAplot"))) {
  dir.create(paste0("./output/", ver, "/MAplot"))
  print(paste("Output MAplot directory for version", ver, "created"))
} else {
  print("DIRECTORY ALREADY EXISTS!!!")
  print("Further analysis will OVERWRITE existing files")
}
# Run DESeq2 and create MAplot for a comparison
source("./R_script/DESeq2/Imm-PDL1_vs_Imm-CTLA4.R")
source("./R_script/DESeq2/Imm-PDL1_vs_Imm-CTR.R")
source("./R_script/DESeq2/Imm-CTLA4_vs_Imm-CTR.R")
source("./R_script/DESeq2/Is-PDL1_vs_Is-CTLA4.R")


# Experiments -------------------------------------------------------------
# Directories for multi, immune, and islet are created in each script
# Run DESeq2 and create MAplot for a comparison
source("./R_script/experiments/multi.R")
source("./R_script/experiments/immune.R")
source("./R_script/experiments/islet.R")


# 3D Plot -----------------------------------------------------------------
# Scripts save outputs in immune, islet created before
# RUN MANUALLY!!


# Volcano Plot ------------------------------------------------------------
# Create output directory
if (!dir.exists(paste0("./output/", ver, "/volcano"))) {
  dir.create(paste0("./output/", ver, "/volcano"))
  print(paste("Output volcano directory for version", ver, "created"))
} else {
  print("DIRECTORY ALREADY EXISTS!!!")
  print("Further analysis will OVERWRITE existing files")
}
# Run scripts that create volcano plots
source("./R_script/volcano/vol_ImPC.R")
source("./R_script/volcano/vol_ImPctr.R")
source("./R_script/volcano/vol_ImCctr.R")
source("./R_script/volcano/vol_IsPC.R")


# Shrinked Volcano Plot ---------------------------------------------------
# Create output directory
if (!dir.exists(paste0("./output/", ver, "/shrkVolcano"))) {
  dir.create(paste0("./output/", ver, "/shrkVolcano"))
  print(paste("Output shrkVolcano directory for version", ver, "created"))
} else {
  print("DIRECTORY ALREADY EXISTS!!!")
  print("Further analysis will OVERWRITE existing files")
}
# Run scripts that create volcano plots (shrinked)
#source("./R_script/shrk_volcano/shrkVol_ImPC.R")
#source("./R_script/shrk_volcano/shrkVol_ImPctr.R")
#source("./R_script/shrk_volcano/shrkVol_ImCctr.R")
#source("./R_script/shrk_volcano/shrkVol_IsPC.R")


# Annotation --------------------------------------------------------------
# Create output directory
if (!dir.exists(paste0("./output/", ver, "/annotation"))) {
  dir.create(paste0("./output/", ver, "/annotation"))
  print(paste("Output annotation directory for version", ver, "created"))
} else {
  print("DIRECTORY ALREADY EXISTS!!!")
  print("Further analysis will OVERWRITE existing files")
}
# Run scripts that create annotated.csv of gene expressions
source("./R_script/annotation/annot_ImPC.R")
source("./R_script/annotation/annot_ImPctr.R")
source("./R_script/annotation/annot_ImCctr.R")
source("./R_script/annotation/annot_IsPC.R")


# Venn Diagram ------------------------------------------------------------
# Create output directory
if (!dir.exists(paste0("./output/", ver, "/venn"))) {
  dir.create(paste0("./output/", ver, "/venn"))
  print(paste("Output venn directory for version", ver, "created"))
} else {
  print("DIRECTORY ALREADY EXISTS!!!")
  print("Further analysis will OVERWRITE existing files")
}
# Run scripts that create overlap in Immune cells
source("./R_script/venn/venn.R")


# plot with DEG ------------------------------------------------------------
# Scripts save outputs in islet directory created before
# Run scripts that create PCA plots with given DEGs
source("./R_script/DEG_plot/DEG_plot_Im.R")
source("./R_script/DEG_plot/DEG_plot_Is.R")


# 3D Plot with DEGs -------------------------------------------------------
# Scripts save outputs in multi, immune, islet created before
# RUN MANUALLY!!


# GO ----------------------------------------------------------------------
# Create output directory
if (!dir.exists(paste0("./output/", ver, "/GO"))) {
  dir.create(paste0("./output/", ver, "/GO"))
  print(paste("Output GO directory for version", ver, "created"))
} else {
  print("DIRECTORY ALREADY EXISTS!!!")
  print("Further analysis will OVERWRITE existing files")
}
# Run scripts that create GO analysis results
source("./R_script/GO/GO_ImPC.R")
# Scripts below are runnable when modified a little, but not required
#source("./R_script/GO/GO_ImPctr.R")
#source("./R_script/GO/GO_ImCctr.R")
#source("./R_script/GO/GO_IsPC.R")


# GO Barplot --------------------------------------------------------------
# Create output directory
if (!dir.exists(paste0("./output/", ver, "/GO_barplot"))) {
  dir.create(paste0("./output/", ver, "/GO_barplot"))
  print(paste("Output GO_barplot directory for version", ver, "created"))
} else {
  print("DIRECTORY ALREADY EXISTS!!!")
  print("Further analysis will OVERWRITE existing files")
}
# Run scripts that create GO analysis results
source("./R_script/GO_barplot/goBar_ImPC.R")


# Gene Location -----------------------------------------------------------


# Violin ------------------------------------------------------------------
# Create output directory
if (!dir.exists(paste0("./output/", ver, "/expression"))) {
  dir.create(paste0("./output/", ver, "/expression"))
  print(paste("Output violin directory for version", ver, "created"))
} else {
  print("DIRECTORY ALREADY EXISTS!!!")
  print("Further analysis will OVERWRITE existing files")
}
# Run scripts that create GO analysis results
source("./R_script/expression/expression_ImPC.R")
source("./R_script/expression/expression_IsPC.R")

