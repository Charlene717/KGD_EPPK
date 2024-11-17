##### Presetting ######
rm(list = ls()) # 清除變數
memory.limit(150000)

#### Load Packages ####
if(!require('Seurat')) {install.packages('Seurat'); library(Seurat)}
if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}
if(!require('stringr')) {install.packages('stringr'); library(stringr)}
if(!require('R.utils')) {install.packages('R.utils'); library(R.utils)}

## 設定包含 .tar 檔案的目錄路徑
# data_dir <- "C:/Users/q2330/Dropbox/KGD_Lab/20241029 (待處理)(EPPK)_線上正常皮膚scRNA-seq/GSE202352_RAW"
data_dir <- "D:/Dropbox/KGD_Lab/20241029 (待處理)(EPPK)_線上正常皮膚scRNA-seq/GSE202352_RAW"

# 列出目錄中的所有 .tar 檔案
tar_files <- list.files(data_dir, pattern = "\\.tar$|\\.tar\\.gz$", full.names = TRUE)

# 確認檔案是否正確列出
if(length(tar_files) == 0) {
  stop("No .tar or .tar.gz files found. Please check the directory path or file extension.")
}

# 提取檔名的 metadata
extract_metadata <- function(filename) {
  base_name <- basename(filename)
  base_name <- sub(".filtered_feature_bc_matrix.tar.gz|.filtered_feature_bc_matrix.tar", "", base_name)
  parts <- unlist(str_split(base_name, "_"))
  
  metadata <- list(
    FullName = base_name,
    GSM_ID = parts[1],
    Sample_Subject = paste(parts[2], parts[3], parts[4], sep = "_"),
    Tissue = parts[5]
  )
  
  return(metadata)
}

# 初始化 Seurat 物件列表
seurat_list <- list()

# 迴圈處理每個檔案，並建立 Seurat 物件
for (file in tar_files) {
  metadata <- extract_metadata(file)
  extract_dir <- file.path(data_dir, metadata$FullName)
  dir.create(extract_dir, showWarnings = FALSE)
  untar(file, exdir = extract_dir)
  
  all_files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
  tenx_files <- list(matrix = NULL, features = NULL, barcodes = NULL)
  
  for (f in all_files) {
    if (grepl("matrix.mtx.gz$", f) && !file.exists(file.path(extract_dir, "matrix.mtx.gz"))) {
      file.copy(f, file.path(extract_dir, "matrix.mtx.gz"), overwrite = TRUE)
      tenx_files$matrix <- file.path(extract_dir, "matrix.mtx.gz")
    } else if (grepl("features.tsv.gz$", f) && !file.exists(file.path(extract_dir, "features.tsv.gz"))) {
      file.copy(f, file.path(extract_dir, "features.tsv.gz"), overwrite = TRUE)
      tenx_files$features <- file.path(extract_dir, "features.tsv.gz")
    } else if (grepl("barcodes.tsv.gz$", f) && !file.exists(file.path(extract_dir, "barcodes.tsv.gz"))) {
      file.copy(f, file.path(extract_dir, "barcodes.tsv.gz"), overwrite = TRUE)
      tenx_files$barcodes <- file.path(extract_dir, "barcodes.tsv.gz")
    }
  }
  
  if (!is.null(tenx_files$matrix) && !is.null(tenx_files$features) && !is.null(tenx_files$barcodes)) {
    message(paste("Loading data from directory:", extract_dir))
    
    tryCatch({
      seurat_data <- Read10X(data.dir = extract_dir)
      seurat_obj <- CreateSeuratObject(counts = seurat_data, project = metadata$FullName)
      
      # 為每個 Seurat 物件進行 NormalizeData 和 FindVariableFeatures
      seurat_obj <- NormalizeData(seurat_obj)
      seurat_obj <- FindVariableFeatures(seurat_obj)
      
      # 添加 metadata
      for (meta_name in names(metadata)) {
        seurat_obj[[meta_name]] <- metadata[[meta_name]]
      }
      
      # 確保細胞名稱唯一化
      seurat_obj <- RenameCells(seurat_obj, add.cell.id = metadata$GSM_ID)
      
      # 將 Seurat 物件儲存在列表中
      seurat_list[[metadata$FullName]] <- seurat_obj
      
    }, error = function(e) {
      message(paste("Error while loading data from:", extract_dir))
      message("Error details:", e)
    })
    
  } else {
    warning(paste("Required compressed files not found for file:", file))
    next
  }
}

# 批次效應整合
if(length(seurat_list) > 0) {
  # 使用 FindIntegrationAnchors 和 IntegrateData 來進行批次效應處理
  features <- SelectIntegrationFeatures(object.list = seurat_list)
  anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
  combined_seurat <- IntegrateData(anchorset = anchors)
  
  # 正規化和降維
  DefaultAssay(combined_seurat) <- "integrated"
  combined_seurat <- ScaleData(combined_seurat)
  combined_seurat <- RunPCA(combined_seurat)
  combined_seurat <- FindNeighbors(combined_seurat, dims = 1:10)
  combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)
  combined_seurat <- RunUMAP(combined_seurat, dims = 1:10)
  
  # 更新 Tissue 資訊
  combined_seurat$Tissue <- ifelse(
    combined_seurat$GSM_ID %in% c("GSM6111848", "GSM6111852"), "hip_dermis",
    ifelse(
      combined_seurat$GSM_ID %in% c("GSM6111849", "GSM6111853"), "palm_dermis",
      ifelse(
        combined_seurat$GSM_ID %in% c("GSM6111850", "GSM6111854"), "hip_epi",
        ifelse(combined_seurat$GSM_ID %in% c("GSM6111851", "GSM6111855"), "palm_epi", combined_seurat$Tissue)
      )
    )
  )
  
  # 保存整合後的 Seurat 物件
  saveRDS(combined_seurat, file = "combined_seurat_integrated.rds")
  
} else {
  stop("No valid Seurat objects were created. Please check file structure and contents.")
}

# 後續分析
DefaultAssay(combined_seurat) <- "RNA"

Marker.lt <- c("EPGN", "EGFR", "KRT1", "KRT9", 
               "KIT", "EP300", "NF1", "COL1A1")
VlnPlot(combined_seurat, features = Marker.lt, ncol = 4)
VlnPlot(combined_seurat, features = Marker.lt[1:4], ncol = 2, group.by = "blueprint.main")
VlnPlot(combined_seurat, features = Marker.lt[5:8], ncol = 2, group.by = "blueprint.main")
FeaturePlot(combined_seurat, features = Marker.lt, ncol = 4)

DimPlot(combined_seurat, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "GSM_ID") + NoLegend()
DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "GSM_ID")
DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "Tissue" ) 
DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "seurat_clusters" ) 
