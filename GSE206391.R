##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

# 檢查並載入必要的套件
if(!require('Seurat')) {install.packages('Seurat'); library(Seurat)}
if(!require('hdf5r')) {install.packages('hdf5r'); library(hdf5r)}
if(!require('Matrix')) {install.packages('Matrix'); library(Matrix)}

# 指定檔案路徑
file_path <- "C:/Charlene/Dataset_Online/GSE206391/GSE206391_Preprocessed_data.h5"

# 打開 HDF5 檔案
h5_file <- H5File$new(file_path, mode = "r")

# 提取表達矩陣數據（在 /layers/counts）
counts_dataset <- h5_file[["layers/counts"]]
dims <- c(16685, 59319)  # 根據檢查結果手動設定維度

# 使用索引提取數據
expression_data <- counts_dataset[1:dims[1], 1:dims[2]]

# 提取基因和細胞名稱
gene_names <- h5_file[["var/gene_name"]][]
cell_names <- h5_file[["obs/_index"]][]

# 確保數據的維度名稱正確對應
rownames(expression_data) <- gene_names
colnames(expression_data) <- cell_names

# 將數據轉換為稀疏矩陣
sparse_matrix <- Matrix::Matrix(expression_data, sparse = TRUE)

# 提取元數據：遍歷 obs 群組中的所有數據集
# 查看 obs 群組內的結構
obs_group <- h5_file[["obs"]]
print(obs_group)  # 確認群組結構

# 列出 obs 群組中的所有對象
dataset_names <- names(obs_group)
print(dataset_names)  # 確認有哪些數據集







# 提取 obs 群組的數據
obs_group <- h5_file[["obs"]]
metadata <- list()

for (dataset_name in names(obs_group)) {
  dataset <- obs_group[[dataset_name]]
  if (inherits(dataset, "H5D")) {  # 確認是否為數據集
    metadata[[dataset_name]] <- dataset[]
  }
}

# 將數據轉為 data.frame
metadata <- as.data.frame(metadata)
print(head(metadata))

# 檢查 __categories 群組（如果存在）
if ("__categories" %in% names(obs_group)) {
  categories_group <- obs_group[["__categories"]]
  for (category_name in names(categories_group)) {
    if (category_name %in% colnames(metadata)) {
      # 提取類別映射表
      category_mapping <- categories_group[[category_name]][]
      
      # 替換數字為字符標籤
      metadata[[category_name]] <- factor(metadata[[category_name]], levels = seq_along(category_mapping), labels = category_mapping)
    }
  }
}

metadata <- as.data.frame(metadata)
head(metadata) # 查看處理後的元數據

# 創建 Seurat 物件
seurat_object <- CreateSeuratObject(counts = sparse_matrix, project = "GSE206391", meta.data = metadata)

# 查看更新後的 Seurat 物件
print(seurat_object)

# 進行後續分析
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
seurat_object <- RunUMAP(seurat_object, dims = 1:10)

# 視覺化結果
DimPlot(seurat_object, reduction = "umap")
DimPlot(seurat_object, reduction = "umap", group.by = "patient")






# # Visualize UMAP results
# DimPlot(acral_melanoma_integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# DimPlot(acral_melanoma_integrated, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "curated.cell.types") + NoLegend()
# DimPlot(acral_melanoma_integrated, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "sample") + NoLegend()
# 
# # Find markers for each cluster
# cluster_markers <- FindAllMarkers(acral_melanoma_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# # Display marker genes
# print(cluster_markers)
# 
# DefaultAssay(acral_melanoma_integrated) <- "RNA"
# 
# Marker.lt <- c("EPGN", "EGFR","KRT1","KRT16",
#                "BRAF","NRAS","NF1","TERT")
# 
# Marker.lt <- c("EPGN", "EGFR","KRT1","KRT9",
#                "KIT","EP300","NF1","COL1A1")
# VlnPlot(acral_melanoma_integrated, features = Marker.lt,ncol = 4)
# VlnPlot(acral_melanoma_integrated, features = Marker.lt[1:4],ncol = 2, group.by = "blueprint.main")
# VlnPlot(acral_melanoma_integrated, features = Marker.lt[5:8],ncol = 2, group.by = "blueprint.main")
# FeaturePlot(acral_melanoma_integrated, features = Marker.lt,ncol = 4)
# 
# FeaturePlot(acral_melanoma_integrated, features = c("EPGN","EGFR","NPEPL1","SERPINA12"),ncol = 4)
# 
# # Save Seurat object
# saveRDS(acral_melanoma_integrated, file = paste0(folder_path, "GSE215121_Acral_Melanoma_SeuratObject_Integrated.rds"))
