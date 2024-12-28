##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
if(!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if(!require('Seurat')) install.packages('Seurat'); library(Seurat)
if(!require('hdf5r')) install.packages("hdf5r"); library(hdf5r)
if(!require('dplyr')) install.packages("dplyr"); library(dplyr)


#### Load Data ####
# folder_path <- "D:/Dropbox/KGD_Lab/20241015 (待整理)(線上資料整理)_Acral Melanoma/Acral melanoma datasets/"
folder_path <- "C:/Users/q2330/Dropbox/KGD_Lab/20241015 (待整理)(線上資料整理)_Acral Melanoma/Acral melanoma datasets/"

file_path <- paste0(folder_path,"GSE189889_acral_2101.debatched.iscva.h5")
  
data <- H5File$new(file_path, mode = "r")

# 列出數據集的結構
data_structure <- data$ls()
print(data_structure)

# 確定矩陣的位置並讀取表達矩陣
# 根據數據結構，更新矩陣的讀取方式，例如：如果矩陣存儲在 "matrix/data"
data_matrix_group <- data[["matrix"]]

# 列出 matrix 群組的結構，以確定數據的位置
matrix_structure <- data_matrix_group$ls()
print(matrix_structure)

# 根據結構讀取表達矩陣 (假設數據在 "data" 鍵中)
# 讀取矩陣形狀信息
matrix_shape <- data_matrix_group[["shape"]]$read()

# 列出 matrix 群組的結構，以確定數據的位置
matrix_structure <- data_matrix_group$ls()
print(matrix_structure)

library(Matrix)

# 讀取稀疏矩陣的數據
indices <- data_matrix_group[["indices"]]$read()
indptr <- data_matrix_group[["indptr"]]$read()
data_values <- data_matrix_group[["data"]]$read()

# 重建稀疏矩陣
expr_matrix <- sparseMatrix(i = indices + 1, p = indptr, x = data_values, dims = c(matrix_shape[1], matrix_shape[2]))

# 確保矩陣有行名和列名
gene_names <- data_matrix_group[["gene_names"]]$read()
cell_names <- data_matrix_group[["barcodes"]]$read()

# 添加行名和列名
if (!is.null(dim(expr_matrix)) && dim(expr_matrix)[1] == length(gene_names) && dim(expr_matrix)[2] == length(cell_names)) {
  rownames(expr_matrix) <- gene_names
  colnames(expr_matrix) <- cell_names
} else {
  stop("基因名稱或細胞名稱與表達矩陣的維度不匹配，請檢查數據讀取過程。")
}


#### Meta data ####
# 列出 artifacts 群組的結構，以確定是否有 meta.data
artifacts_group <- data[["artifacts"]]
artifacts_structure <- artifacts_group$ls()
print(artifacts_structure)


# 列出 artifacts/all 群組的結構
all_group <- artifacts_group[["all"]]
all_structure <- all_group$ls()
print(all_structure)

# # 列出 artifacts/lymphoids 群組的結構
# lymphoids_group <- artifacts_group[["lymphoids"]]
# lymphoids_structure <- lymphoids_group$ls()
# print(lymphoids_structure)
# 
# # 列出 artifacts/myeloids 群組的結構
# myeloids_group <- artifacts_group[["myeloids"]]
# myeloids_structure <- myeloids_group$ls()
# print(myeloids_structure)



# 讀取 artifacts/all 中的 meta.data
covs <- all_group[["covs"]]$read()
pcs <- all_group[["pcs"]]$read()

# 將 covs 和 pcs 合併到 meta.data
meta_data <- as.data.frame(cbind(covs, pcs))
rownames(meta_data) <- cell_names

# 創建 Seurat 對象並添加 meta.data
acral_melanoma <- CreateSeuratObject(counts = expr_matrix, project = "AcralMelanoma", min.cells = 3, min.features = 200, meta.data = meta_data)




# 添加質控指標
acral_melanoma[["percent.mt"]] <- PercentageFeatureSet(acral_melanoma, pattern = "^MT-")

# 質控過濾
acral_melanoma <- subset(acral_melanoma, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# 標準化數據
acral_melanoma <- NormalizeData(acral_melanoma)

# 識別高變基因
acral_melanoma <- FindVariableFeatures(acral_melanoma, selection.method = "vst", nfeatures = 2000)

# 整合數據集
all.genes <- rownames(acral_melanoma)
acral_melanoma <- ScaleData(acral_melanoma, features = all.genes)

# 主成分分析 (PCA)
acral_melanoma <- RunPCA(acral_melanoma, features = VariableFeatures(object = acral_melanoma))

# 視覺化 PCA 結果
ElbowPlot(acral_melanoma)

# 使用 PCA 進行聚類分析
acral_melanoma <- FindNeighbors(acral_melanoma, dims = 1:10)
acral_melanoma <- FindClusters(acral_melanoma, resolution = 0.5)

# UMAP 降維
acral_melanoma <- RunUMAP(acral_melanoma, dims = 1:10)

# 視覺化 UMAP 結果
DimPlot(acral_melanoma, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(acral_melanoma, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "curated.cell.types" ) + NoLegend()
DimPlot(acral_melanoma, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "blueprint.main" ) 
DimPlot(acral_melanoma, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "sample" ) + NoLegend()

# 找出每個群集的標誌基因
cluster_markers <- FindAllMarkers(acral_melanoma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(acral_melanoma, features = top10$gene) + NoLegend()


# 顯示標誌基因
print(cluster_markers)

VlnPlot(acral_melanoma, features = c("EPGN", "EGFR"))


Marker.lt <- c("EPGN", "EGFR","KRT1","KRT16",
               "BRAF","NRAS","NF1","TERT")

VlnPlot(acral_melanoma, features = Marker.lt,ncol = 4)
FeaturePlot(acral_melanoma, features = Marker.lt,ncol = 4)

VlnPlot(acral_melanoma, features = Marker.lt, ncol = 4, group.by = "blueprint.main") 

VlnPlot(acral_melanoma, features = Marker.lt[1:4], ncol = 2, group.by = "blueprint.main") 
VlnPlot(acral_melanoma, features = Marker.lt[5:8], ncol = 2, group.by = "blueprint.main")
VlnPlot(acral_melanoma, features = Marker.lt[5:8], ncol = 2, group.by = "blueprint.main", pt.size = 0, alpha = 0.3)
VlnPlot(acral_melanoma, features = Marker.lt[5:8], ncol = 2, group.by = "blueprint.main", pt.size = 0.01, alpha = 0.1)
acral_melanoma$blueprint.main %>% table()

# VlnPlot(acral_melanoma, features = Marker.lt, ncol = 2, group.by = "blueprint.main") 

# library(ggplot2)
# 
# plots <- VlnPlot(acral_melanoma, features = Marker.lt, ncol = 4, group.by = "blueprint.main", combine = FALSE)
# plots <- lapply(plots, function(plot) {
#   plot + theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))
# })
# CombinePlots(plots = plots, ncol = 4)



# 保存 Seurat 對象
saveRDS(acral_melanoma, file = paste0(folder_path, "GSE189889_Acral_Melanoma_SeuratObject.rds"))


#################################################################################
acral_melanoma <- readRDS("C:/Users/q2330/Dropbox/KGD_Lab/20241015 (待整理)(線上資料整理)_Acral Melanoma/Acral melanoma datasets/GSE189889_Acral_Melanoma_SeuratObject.rds")
VlnPlot(acral_melanoma, features = c("EPGN", "EGFR"))