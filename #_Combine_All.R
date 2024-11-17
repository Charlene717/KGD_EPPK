##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
if(!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if(!require('Seurat')) install.packages('Seurat'); library(Seurat)
if(!require('hdf5r')) install.packages("hdf5r"); library(hdf5r)
if(!require('dplyr')) install.packages("dplyr"); library(dplyr)
if(!require('Matrix')) install.packages("Matrix"); library(Matrix)

Marker.lt <- c("EPGN", "EGFR", "KRT1", "KRT9", 
               "KIT", "EP300", "NF1", "COL1A1")

#### GSE215121 ####
folder_path <- "C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/"
# load("C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/GSE215121_Acral_Melanoma_Temp.RData")
# DimPlot(acral_melanoma, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# 
# # Save Seurat object
# saveRDS(acral_melanoma, file = paste0(folder_path, "GSE215121_Acral_Melanoma_SeuratObject.rds"))

GSE215121_Seurat <-  readRDS(file ="C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/GSE215121_Acral_Melanoma_SeuratObject.rds")
DimPlot(GSE215121_Seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

GSE215121_Seurat$Dataset <- "GSE215121"
DimPlot(GSE215121_Seurat, reduction = "umap",group.by = "Dataset", label = TRUE, pt.size = 0.5) + NoLegend()

VlnPlot(GSE215121_Seurat, features = Marker.lt, ncol = 4)


#### KGD_EPPK ####
# load("C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/#_KGL_Lab_EPPK_HS23378_2862.Temp.RData")
# # Save Seurat object
# saveRDS(seurat_obj, file = paste0(folder_path, "KGL_Lab_EPPK_HS23378_2862.rds"))

KGD_EPPK__HS23378_2862_Seurat <-  readRDS(file ="C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/KGL_Lab_EPPK_HS23378_2862.rds")
DimPlot(KGD_EPPK__HS23378_2862_Seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
KGD_EPPK__HS23378_2862_Seurat$Dataset <- "KGD_EPPK__HS23378_2862"
DimPlot(KGD_EPPK__HS23378_2862_Seurat, reduction = "umap",group.by = "Dataset", label = TRUE, pt.size = 0.5) + NoLegend()


VlnPlot(KGD_EPPK__HS23378_2862_Seurat, features = Marker.lt, ncol = 4)



