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


seurat_object$patient %>% table()
seurat_object$replicates %>% table()


#### Add disease ####
# 定義各疾病對應的 patient 編號
ad_patients <- c(2, 5, 8, 11, 15, 20, 34, 35, 36)
lp_patients <- c(3, 6, 9, 12, 14, 25, 26, 27, 28, 30, 37)
pso_patients <- c(1, 10, 13, 19, 22, 29, 31, 32, 33)

# 新增 disease 欄位
seurat_object$disease <- NA  # 初始化為 NA

# 設定規則
seurat_object$disease[seurat_object$patient %in% ad_patients] <- "AD"
seurat_object$disease[seurat_object$patient %in% lp_patients] <- "LP"
seurat_object$disease[seurat_object$patient %in% pso_patients] <- "Pso"

# 確認新增的欄位
table(seurat_object$disease, useNA = "ifany")  # 檢查各類別的分佈

DimPlot(seurat_object, reduction = "umap", group.by = "disease")

seurat_object$patient %>% table()
seurat_object$patient %>% unique()


# 新增 group 欄位
seurat_object$group <- paste0(seurat_object$disease, "_", sub(".*_", "", seurat_object$replicates))
seurat_object$groupLN <- paste0(sub(".*_", "", seurat_object$replicates))

# 確認結果
table(seurat_object$group, useNA = "ifany")  # 查看每個 group 的分佈
DimPlot(seurat_object, reduction = "umap", group.by = "group")


#### Remove NA ####
# 提取細胞元數據
metadata <- seurat_object@meta.data

# 找到 patient 為 NA 的條目
valid_cells <- !is.na(metadata$patient)

# 根據有效細胞更新 Seurat 物件
seurat_object <- seurat_object[, valid_cells]

# 確認結果
seurat_object$patient %>% unique()  # 查看是否還有 NA

DimPlot(seurat_object, reduction = "umap", group.by = "group")

seurat_object$group %>% table()

Marker.lt <- c("EPGN", "EGFR","KRT1","KRT9",
               "KIT","EP300","NF1","COL1A1")
VlnPlot(seurat_object, features = Marker.lt[1:4],ncol = 2, group.by = "group")





# 檢查並載入必要套件
if (!require('ggplot2')) { install.packages('ggplot2'); library(ggplot2) }
if (!require('reshape2')) { install.packages('reshape2'); library(reshape2) }

# 定義要繪製的 Marker
Marker.lt <- c("EPGN", "EGFR")

# 提取數據
plot_data <- FetchData(seurat_object, vars = c(Marker.lt, "group", "groupLN"))

# 將數據轉換為長格式以適應 ggplot2
plot_data_long <- reshape2::melt(plot_data, id.vars = c("group", "groupLN"), variable.name = "Marker", value.name = "Expression")

# 繪製 Boxplot
ggplot(plot_data_long, aes(x = group, y = Expression, fill = groupLN)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +  # 繪製透明度更高的箱線圖
  geom_jitter(width = 0.2, size = 1, alpha = 0.2, shape = 21, stroke = 0.5, color = "black") +  # 添加透明點
  facet_wrap(~ Marker, ncol = 2, scales = "free_y") +  # 分面顯示每個 Marker
  labs(x = "Group", y = "Expression Level", fill = "GroupLN") +  # 添加標籤
  theme_minimal(base_size = 15) +  # 使用簡潔主題並設定基礎文字大小
  theme(
    strip.text = element_text(size = 24, face = "bold"),  # 分面標題加大
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # X 軸標籤加大並旋轉
    axis.text.y = element_text(size = 22),  # Y 軸標籤加大
    axis.title.x = element_text(size = 18, face = "bold"),  # X 軸標題加大
    axis.title.y = element_text(size = 18, face = "bold"),  # Y 軸標題加大
    legend.title = element_text(size = 18, face = "bold"),  # 圖例標題加大
    legend.text = element_text(size = 16)  # 圖例文字加大
  )



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
