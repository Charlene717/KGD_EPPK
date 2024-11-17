##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
if(!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if(!require('Seurat')) install.packages('Seurat'); library(Seurat)
if(!require('hdf5r')) install.packages("hdf5r"); library(hdf5r)
if(!require('dplyr')) install.packages("dplyr"); library(dplyr)
if(!require('Matrix')) install.packages("Matrix"); library(Matrix)


#### GSE215121 ####
folder_path <- "C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/"
load("C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/GSE215121_Acral_Melanoma_Temp.RData")
DimPlot(acral_melanoma, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Save Seurat object
saveRDS(acral_melanoma, file = paste0(folder_path, "GSE215121_Acral_Melanoma_SeuratObject.rds"))

GSE215121_Seurat <-  readRDS(file ="C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/GSE215121_Acral_Melanoma_SeuratObject.rds")
DimPlot(GSE215121_Seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



