##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
if(!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if(!require('Seurat')) install.packages('Seurat'); library(Seurat)
if(!require('hdf5r')) install.packages("hdf5r"); library(hdf5r)
if(!require('dplyr')) install.packages("dplyr"); library(dplyr)

#### Load data ####

Set_Dir <- "F:/20230522 fixed single cells 10X_2795+2818+2862+2901/20230522 fixed single cells 10X_2795+2818+2862+2901/results_HS23378_2862/per_sample_outs/2862/count/sample_filtered_feature_bc_matrix"
seurat_data <- Read10X(data.dir = Set_Dir)
seurat_obj <- CreateSeuratObject(counts = seurat_data)


# 添加質控指標
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# 質控過濾
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# 標準化數據
seurat_obj <- NormalizeData(seurat_obj)

# 識別高變基因
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 整合數據集
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

# 主成分分析 (PCA)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# 視覺化 PCA 結果
ElbowPlot(seurat_obj)

# 使用 PCA 進行聚類分析
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# UMAP 降維
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# 視覺化 UMAP 結果
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)

DefaultAssay(seurat_obj) <- "RNA"


Marker.lt <- c("EPGN", "EGFR","KRT1","KRT9",
               "KIT","EP300","NF1","COL1A1")
VlnPlot(seurat_obj, features = Marker.lt,ncol = 4)
VlnPlot(seurat_obj, features = Marker.lt[1:4],ncol = 2, group.by = "blueprint.main")
VlnPlot(seurat_obj, features = Marker.lt[5:8],ncol = 2, group.by = "blueprint.main")
FeaturePlot(seurat_obj, features = Marker.lt,ncol = 4)

FeaturePlot(seurat_obj, features = c("EPGN","EGFR","NPEPL1","SERPINA12"),ncol = 4)
