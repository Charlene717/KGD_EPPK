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



#### KGD_EPPK ####
# load("C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/#_KGL_Lab_EPPK_HS23378_2862.Temp.RData")
# # Save Seurat object
# saveRDS(seurat_obj, file = paste0(folder_path, "KGL_Lab_EPPK_HS23378_2862.rds"))

KGD_EPPK_HS23378_2862_Seurat <-  readRDS(file ="C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/KGL_Lab_EPPK_HS23378_2862.rds")
DimPlot(KGD_EPPK_HS23378_2862_Seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
KGD_EPPK_HS23378_2862_Seurat$Dataset <- "KGD_EPPK_HS23378_2862"
DimPlot(KGD_EPPK_HS23378_2862_Seurat, reduction = "umap",group.by = "Dataset", label = TRUE, pt.size = 0.5) + NoLegend()


VlnPlot(KGD_EPPK_HS23378_2862_Seurat, features = Marker.lt, ncol = 4)


#### GSE215121 Acral melanoma ####
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


#### GSE189889 Acral melanoma ####
GSE189889_Seurat <-  readRDS(file ="C:/Users/q2330/Dropbox/KGD_Lab/20241015 (已整理)(EPPK)(線上資料整理)_Acral Melanoma/Acral melanoma datasets/GSE189889_Acral_Melanoma_SeuratObject.rds")
DimPlot(GSE189889_Seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
GSE189889_Seurat$Dataset <- "GSE189889"
DimPlot(GSE189889_Seurat, reduction = "umap",group.by = "Dataset", label = TRUE, pt.size = 0.5) + NoLegend()


VlnPlot(GSE189889_Seurat, features = Marker.lt, ncol = 4)

#### GSE202352 Normal  ####
GSE202352_Normal <-  readRDS(file ="C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/GSE202352_combined_seurat_integrated.rds")
DimPlot(GSE202352_Normal, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
GSE202352_Normal$Dataset <- "GSE202352"
DimPlot(GSE202352_Normal, reduction = "umap",group.by = "Dataset", label = TRUE, pt.size = 0.5) + NoLegend()


VlnPlot(GSE202352_Normal, features = Marker.lt, ncol = 4)


###############################################################################

DefaultAssay(KGD_EPPK_HS23378_2862_Seurat) <- "RNA"
DefaultAssay(GSE215121_Seurat) <- "RNA"
DefaultAssay(GSE189889_Seurat) <- "RNA"
DefaultAssay(GSE202352_Normal) <- "RNA"

#### Preprocessing: Find Integration Anchors ####
# List of Seurat objects
seurat.list <- list(KGD_EPPK_HS23378_2862_Seurat, GSE215121_Seurat, GSE189889_Seurat, GSE202352_Normal)

# Standardize the preprocessing steps
for (i in 1:length(seurat.list)) {
  seurat.list[[i]] <- NormalizeData(seurat.list[[i]])
  seurat.list[[i]] <- FindVariableFeatures(seurat.list[[i]], selection.method = "vst", nfeatures = 2000)
  seurat.list[[i]] <- ScaleData(seurat.list[[i]], features = VariableFeatures(object = seurat.list[[i]]), verbose = FALSE)
}

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:30)

#### Integrate Data ####
integrated.data <- IntegrateData(anchorset = anchors, dims = 1:30)

# Switch to integrated assay for downstream analysis
DefaultAssay(integrated.data) <- "integrated"

#### Run Standard Workflow ####
# Scaling, PCA, UMAP, and Clustering
integrated.data <- ScaleData(integrated.data, features = rownames(integrated.data), verbose = FALSE)
integrated.data <- RunPCA(integrated.data, npcs = 30, verbose = FALSE)
integrated.data <- RunUMAP(integrated.data, reduction = "pca", dims = 1:30)
integrated.data <- FindNeighbors(integrated.data, dims = 1:30)
integrated.data <- FindClusters(integrated.data, resolution = 0.5)

#### Visualization ####
# UMAP Plot
DimPlot(integrated.data, reduction = "umap", group.by = "Dataset", label = TRUE, pt.size = 0.5) + NoLegend()

# Violin Plot for Marker Genes
VlnPlot(integrated.data, features = Marker.lt, ncol = 4)
