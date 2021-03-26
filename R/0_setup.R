#set up folders
dir.create("data")
dir.create("data/human_counts")
dir.create("data/murine_counts")
dir.create("R")
dir.create("plots")
dir.create("plots/heatmaps")
dir.create("plots/umap")
dir.create("plots/others")


#for these analyses you will need the following packages
#install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

#the following code was adjusted from url https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c(
  "Seurat",
  "assertthat",
  "tools",
  "tidyverse",
  "data.table",
  "readr",
  "RaceID",
  "Seurat",
  "Matrix",
  "clustree",
  "cowplot",
  "assertthat",
  "viridis",
  "readxl",
  "RColorBrewer",
  "pheatmap",
  "MASS",
  "rebus",
  "emmeans"
)


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
