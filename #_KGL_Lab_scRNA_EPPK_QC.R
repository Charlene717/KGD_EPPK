## Ref: https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
if(!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if(!require('Seurat')) install.packages('Seurat'); library(Seurat)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('patchwork')) install.packages('patchwork'); library(patchwork)

# if (!requireNamespace("Seurat", quietly = TRUE)) { install.packages("Seurat") }; library(Seurat)

#### Set Parameter ####
Set_QC <- TRUE
# Set export parameters
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 10) # Generate unique time-based ID
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))

Set_note <- paste0(Name_FileID,"_QC",Set_QC)

Name_Export <- paste0(Name_FileID)
Name_ExportFolder <- paste0("Export_",Set_note)
# Create export folder if it does not exist
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}

Name_ExportFolder_QC <- paste0(Name_ExportFolder,"/QC")
# Create export folder if it does not exist
if (!dir.exists(Name_ExportFolder_QC)){dir.create(Name_ExportFolder_QC)}


#### Load Data ####
# 設定資料夾路徑
main_folder <- "C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/KGD_EPPK/Input/10X/"

# 獲取資料夾中的所有子資料夾
sub_folders <- list.dirs(main_folder, full.names = TRUE, recursive = FALSE)

# 初始化空的list來存放各個Seurat物件
seurat_list_before <- list()
seurat_list <- list()

# 逐一讀取每個子資料夾
for (folder in sub_folders) {
  
  # 獲取子資料夾名稱
  folder_name <- basename(folder)
  
  # 根據子資料夾名稱提取檔案編號和組織名稱
  # 這裡假設子資料夾名稱的格式是類似 "TN254_DDEB"
  parts <- unlist(strsplit(folder_name, "_"))
  sample_id <- parts[1]
  tissue_type <- parts[2]
  
  # 讀取10x檔案
  data_dir <- file.path(folder)
  if (file.exists(file.path(data_dir, "barcodes.tsv.gz")) &&
      file.exists(file.path(data_dir, "features.tsv.gz")) &&
      file.exists(file.path(data_dir, "matrix.mtx.gz"))) {
    
    # 使用Read10X來讀取資料
    counts <- Read10X(data.dir = data_dir)
    
    # 如果 counts 是列表，選取 Gene Expression 數據
    if (is.list(counts)) {
      counts <- counts[["Gene Expression"]]
    }
    
    # 建立Seurat物件
    seurat_obj <- CreateSeuratObject(counts = counts, project = folder_name)
    
    # 增加Metadata
    seurat_obj$SampleID <- sample_id
    seurat_obj$TissueType <- tissue_type
    
    # 計算線粒體基因百分比
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
    
    # 計算基因複雜性
    seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
    seurat_obj@meta.data$log10GenesPerUMI <- seurat_obj$log10GenesPerUMI
    
    # 保存QC前的Seurat物件作為比較基礎
    seurat_obj_before <- seurat_obj
    seurat_list_before[[folder_name]] <- seurat_obj_before
    
    # 繪製過濾標準的分布圖來檢查數據分布
    plot1 <- VlnPlot(seurat_obj_before, features = "nCount_RNA", pt.size = 0.1) + ggtitle("nCount_RNA Distribution")
    plot2 <- VlnPlot(seurat_obj_before, features = "nFeature_RNA", pt.size = 0.1) + ggtitle("nFeature_RNA Distribution")
    plot3 <- VlnPlot(seurat_obj_before, features = "percent.mt", pt.size = 0.1) + ggtitle("percent.mt Distribution")
    plot4 <- VlnPlot(seurat_obj_before, features = "log10GenesPerUMI", pt.size = 0.1) + ggtitle("log10GenesPerUMI Distribution")
    combined_dist_plot <- (plot1 | plot2) / (plot3 | plot4)
    ggsave(filename = paste0(Name_ExportFolder_QC,"/",folder_name, "_distribution_plots_before_QC.png"), plot = combined_dist_plot, width = 14, height = 10)
    
    # 品質控制過濾
    seurat_obj <- subset(x = seurat_obj,  
                         subset= (nCount_RNA >= 500) & 
                           (nFeature_RNA >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (percent.mt < 20))
    
    # 繪製QC前後的比較圖像（在同一張圖上進行比較）
    # Cell count before and after filtering
    before_count <- nrow(seurat_obj_before@meta.data)
    after_count <- nrow(seurat_obj@meta.data)
    
    p0 <- ggplot() +
      geom_bar(data = data.frame(Condition = rep("Before QC", before_count)), aes(x = Condition, fill = Condition), alpha = 0.3, color = "#5d5e5e", stat = "count") +
      geom_text(aes(x = "Before QC", y = before_count, label = before_count), vjust = -0.5) +  # Add text annotation for "Before QC"
      geom_bar(data = data.frame(Condition = rep("After QC", after_count)), aes(x = Condition, fill = Condition), alpha = 0.3, color = "#01814a", stat = "count") +
      geom_text(aes(x = "After QC", y = after_count, label = after_count), vjust = -0.5) +  # Add text annotation for "After QC"
      ggtitle("Cell count (Before and After QC)") +
      xlab("Condition") +
      ylab("Cell Count") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "white")) +
      scale_fill_manual(values = c("Before QC" = "#5d5e5e", "After QC" = "#01814a")) +
      theme(legend.position = "none") + 
      guides(fill = guide_legend(title = "QC Status"))
    
    # Transcripts per cell before and after filtering
    p1 <- ggplot() +
      geom_density(data = seurat_obj_before@meta.data, aes(x = nCount_RNA, fill = "Before QC"), alpha = 0.3, color = "#5d5e5e") +
      geom_density(data = seurat_obj@meta.data, aes(x = nCount_RNA, fill = "After QC"), alpha = 0.3, color = "#01814a") +
      ggtitle("Transcripts per cell (Before and After QC)") +
      xlab("nUMI") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "white")) +
      scale_fill_manual(values = c("Before QC" = "#5d5e5e", "After QC" = "#01814a"))
    
    # Distribution of genes per cell before and after filtering
    p2 <- ggplot() +
      geom_density(data = seurat_obj_before@meta.data, aes(x = nFeature_RNA, fill = "Before QC"), alpha = 0.3, color = "#5d5e5e") +
      geom_density(data = seurat_obj@meta.data, aes(x = nFeature_RNA, fill = "After QC"), alpha = 0.3, color = "#01814a") +
      ggtitle("Genes per cell (Before and After QC)") +
      xlab("nGene") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "white")) +
      scale_fill_manual(values = c("Before QC" = "#5d5e5e", "After QC" = "#01814a"))
    
    # Mitochondrial counts ratio before and after filtering
    p3 <- ggplot() +
      geom_density(data = seurat_obj_before@meta.data, aes(x = percent.mt, fill = "Before QC"), alpha = 0.3, color = "#5d5e5e") +
      geom_density(data = seurat_obj@meta.data, aes(x = percent.mt, fill = "After QC"), alpha = 0.3, color = "#01814a") +
      ggtitle("Mitochondrial counts ratio (Before and After QC)") +
      xlab("mitoRatio") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "white")) +
      scale_fill_manual(values = c("Before QC" = "#5d5e5e", "After QC" = "#01814a"))
    
    # Complexity before and after filtering
    p4 <- ggplot() +
      geom_density(data = seurat_obj_before@meta.data, aes(x = log10GenesPerUMI, fill = "Before QC"), alpha = 0.3, color = "#5d5e5e") +
      geom_density(data = seurat_obj@meta.data, aes(x = log10GenesPerUMI, fill = "After QC"), alpha = 0.3, color = "#01814a") +
      ggtitle("Complexity (Before and After QC)") +
      xlab("log10GenesPerUMI") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "white")) +
      scale_fill_manual(values = c("Before QC" = "#5d5e5e", "After QC" = "#01814a"))
    
    # 組合所有QC圖像成一張大圖
    combined_plot <- (p0 | p1 | p2 | p3 | p4) + plot_layout(guides = 'collect')
    
    # 保存組合後的QC圖像
    ggsave(filename = paste0(Name_ExportFolder_QC, "/", folder_name, "_QC_comparison_combined.png"), plot = combined_plot, width = 30, height = 6)
    
    # 正規化數據並尋找變異特徵
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj)
    
    # 將Seurat物件加入到列表中
    seurat_list[[folder_name]] <- seurat_obj
  }
}

# # 繪製整合前的 QC 圖像
# if (length(seurat_list_before) > 1) {
#   combined_before <- merge(seurat_list_before[[1]], y = seurat_list_before[-1])
#   combined_before_plot <- VlnPlot(combined_before, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "log10GenesPerUMI"), group.by = "TissueType", pt.size = 0.1, ncol = 2) + ggtitle("Before QC (by Tissue Type)")
#   ggsave(filename = paste0(Name_ExportFolder_QC,"/","combined_before_QC.png"), plot = combined_before_plot, width = 16, height = 8)
# }
# 
# # 繪製整合後的 QC 圖像
# if (length(seurat_list) > 1) {
#   combined_after <- merge(seurat_list[[1]], y = seurat_list[-1])
#   combined_after_plot <- VlnPlot(combined_after, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "log10GenesPerUMI"), group.by = "TissueType", pt.size = 0.1, ncol = 2) + ggtitle("After QC (by Tissue Type)")
#   ggsave(filename = paste0(Name_ExportFolder_QC,"/","combined_after_QC.png"), plot = combined_after_plot, width = 16, height = 8)
# }

# 合併整合前和整合後的數據並重命名細胞名稱以避免重複
if (length(seurat_list_before) > 1 && length(seurat_list) > 1) {
  # 使用 RenameCells 确保唯一性
  seurat_list_before <- lapply(seurat_list_before, function(x) RenameCells(x, add.cell.id = x@project.name))
  seurat_list <- lapply(seurat_list, function(x) RenameCells(x, add.cell.id = x@project.name))
  
  # 合併整合前和整合後的 Seurat 物件
  combined_before <- merge(seurat_list_before[[1]], y = seurat_list_before[-1])
  combined_before$QC_Status <- "Before QC"
  combined_after <- merge(seurat_list[[1]], y = seurat_list[-1])
  combined_after$QC_Status <- "After QC"
  
  # 合併數據
  combined_data <- merge(combined_before, combined_after)
  
  # 手动绘制小提琴图以完全控制颜色和透明度，并正确分离点
  library(ggplot2)
  plot_list <- list()
  features <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "log10GenesPerUMI")
  for (feature in features) {
    p <- ggplot(combined_data@meta.data, aes_string(x = "TissueType", y = feature)) +
      geom_violin(aes(fill = QC_Status), alpha = 0.5, scale = "width", position = position_dodge(width = 0.75)) +
      geom_jitter(aes(color = QC_Status), size = 1, alpha = 0.1, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) + # 使用 position_jitterdodge 并调整大小和透明度
      scale_fill_manual(values = c("Before QC" = "#5d5e5e", "After QC" = "#01814a")) +
      scale_color_manual(values = c("Before QC" = "#5d5e5e", "After QC" = "#01814a")) +
      ggtitle(feature) +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "white"))
    plot_list[[feature]] <- p
  }
  
  # 使用 patchwork 组合图像并将图例放到最下面
  library(patchwork)
  vln_plot_combined <- wrap_plots(plot_list, ncol = 2) + 
    plot_layout(guides = 'collect') & 
    theme(legend.position = "bottom")
  
  # 繪製細胞數量比較的柱狀圖
  cell_count_data <- data.frame(
    TissueType = c(sapply(seurat_list_before, function(x) unique(x$TissueType)), sapply(seurat_list, function(x) unique(x$TissueType))),
    QC_Status = c(rep("Before QC", length(seurat_list_before)), rep("After QC", length(seurat_list))),
    CellCount = c(sapply(seurat_list_before, function(x) nrow(x@meta.data)), sapply(seurat_list, function(x) nrow(x@meta.data)))
  )
  
  # 避免顯示多餘的數字
  cell_count_data <- cell_count_data %>%
    group_by(TissueType, QC_Status) %>%
    summarise(CellCount = sum(CellCount))
  
  cell_count_plot <- ggplot(cell_count_data, aes(x = TissueType, y = CellCount, fill = QC_Status)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = CellCount), position = position_dodge(width = 0.9), vjust = -0.5) +
    ggtitle("Cell Count (Before and After QC) by Tissue Type") +
    xlab("Tissue Type") +
    ylab("Cell Count") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"), 
          plot.background = element_rect(fill = "white", color = NA), 
          legend.position = "bottom") + 
    scale_fill_manual(values = c("Before QC" = "#5d5e5e", "After QC" = "#01814a"))
  
  # 保存 QC 比較圖
  ggsave(filename = paste0(Name_ExportFolder_QC, "/", "combined_QC_comparison_vln.png"), plot = vln_plot_combined, width = 12, height = 10)
  ggsave(filename = paste0(Name_ExportFolder_QC, "/", "cell_count_comparison.png"), plot = cell_count_plot, width = 8, height = 6)
}


#### 整合資料和QC ####
# # 保留 seurat_list，刪除其餘對象
# rm(list = setdiff(ls(), "seurat_list"))
# gc()  # 強制執行垃圾回收以釋放記憶體
# seurat_list[["TN131"]] <- NULL

# 檢查是否有足夠的 Seurat 物件進行整合
if (length(seurat_list) > 1) {
  # 使用 SelectIntegrationFeatures 和 FindIntegrationAnchors 來進行批次效應處理
  features <- SelectIntegrationFeatures(object.list = seurat_list)
  seurat_list <- lapply(seurat_list, function(x) {
    x <- ScaleData(x, features = features)  # 縮放數據
    x <- RunPCA(x, features = features)     # 進行 PCA 分析
    return(x)
  })
  
  # 找到整合的 anchors
  anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
  
  # 整合數據
  combined_seurat <- IntegrateData(anchorset = anchors)
  
  # 正規化整合後的數據
  combined_seurat <- ScaleData(combined_seurat)
  combined_seurat <- RunPCA(combined_seurat)
  combined_seurat <- FindNeighbors(combined_seurat, dims = 1:10)
  combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)
  combined_seurat <- RunUMAP(combined_seurat, dims = 1:10)
  
  ## 保存整合後的 Seurat 物件
  saveRDS(combined_seurat, file = paste0(Name_ExportFolder,"/",Name_Export, "_QC.rds"))
  ## Export RData
  save.image(paste0(Name_ExportFolder,"/", Name_Export,"_QC.RData"))
  
} else {
  stop("No valid Seurat objects were created. Please check file structure and contents.")
}


DimPlot(combined_seurat, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "orig.ident") + NoLegend()
DimPlot(combined_seurat, reduction = "umap", label = F, pt.size = 0.5, group.by = "orig.ident") 
DimPlot(combined_seurat, reduction = "umap", label = F, pt.size = 0.5, group.by = "seurat_clusters") 

ggplot(combined_seurat@meta.data, aes(x = seurat_clusters, fill = orig.ident)) +
  geom_bar(color = "black", size = 1.2) +  # 黑色粗框
  labs(x = "Cell type", y = "Count", fill = "Tissue Type") +
  theme_minimal() +
  scale_x_discrete(drop = FALSE) +  # 確保所有 cluster 類別顯示在 x 軸
  theme(
    axis.text = element_text(size = 14, face = "bold"),       # x, y 軸標籤加大加粗
    axis.title = element_text(size = 16, face = "bold"),      # x, y 軸標題加大加粗
    legend.title = element_text(size = 14, face = "bold"),    # 圖例標題加大加粗
    legend.text = element_text(size = 12, face = "bold")      # 圖例文字加大加粗
  )


combined_seurat <- ScaleData(combined_seurat, do.center = FALSE)  # 避免居中產生負值

DefaultAssay(combined_seurat) <- "RNA"
plots <- VlnPlot(combined_seurat, features = c("EPGN", "EGFR"), split.by = "orig.ident", group.by = "seurat_clusters",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(combined_seurat, features = c("EPGN", "EGFR"), ncol = 4)

################################################################################
#### ROGUE ####
## (T) Define Subcluster
DimPlot(seuratObject, reduction = "umap", group.by = "Cell_Type")
DefaultAssay(seuratObject) <- "integrated"
if(Set_Run_Subcluster){ source("FUN_Subcluster_by_ROGUE.R") }
seuratObject <- Define_Subcluster(seuratObject,
                                  ROGUE_Thr = Set_ROGUE_Thr,
                                  SubClust_Resolution = Set_Clust_resolution,
                                  Num_PCA = Set_Num_PCA)

DimPlot(seuratObject, reduction = "umap", group.by = "Cell_Type")
DimPlot(seuratObject, reduction = "umap", group.by = "seurat_clusters")


################################################################################
#### Cell Type Annotation ####



## Reference
seuratObject_Ref <-  readRDS(file ="C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/KGD_EPPK/Input/rds/Reference/Combined_Fskin_obj_2_1_2_3_webportal_Adult10000.rds")

Set_Ref_Type <-  "Broad"
if(Set_Ref_Type == "Broad"){
  seuratObject_Ref$Actual_Cell_Type <- seuratObject_Ref$Cell_Type_Assigned
}else if(Set_Ref_Type == "Fine"){
  seuratObject_Ref$Actual_Cell_Type <- seuratObject_Ref$annotation_fine_merged
}

DimPlot(seuratObject_Ref, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "Actual_Cell_Type") + NoLegend()
seuratObject_Ref$Actual_Cell_Type %>% as.character() %>% table()
seuratObject_Ref$Actual_Cell_Type %>% is.na() %>% summary()

source("#_FUN_CellTypeAnnot.R")


# 定義標註函數列表
annotation_functions <- list(
  Run_singleR,
  Run_scPred
  # Run_singleR,
  # Run_scmap,
  # Run_SCINA,
  # Run_scPred,
  # Run_CHETAH,
  # Run_scClassify,
  # Run_Seurat_Annot
)



# 依次執行每個標註函數
for (func in annotation_functions) {
  try({ seuratObject_Sample <- func(combined_seurat, seuratObject_Ref) })
}


DimPlot(seuratObject_Sample,group.by = "label_singleR_NoReject", reduction = "umap")
DimPlot(seuratObject_Sample,group.by = "label_singleR", reduction = "umap")


DimPlot(seuratObject_Sample,group.by = "label_scPred_NoReject", reduction = "umap")
DimPlot(seuratObject_Sample,group.by = "label_scPred", reduction = "umap")





################################################################################
# 後續分析
DefaultAssay(combined_seurat) <- "RNA"

Marker.lt <- c("EPGN", "EGFR", "KRT1", "KRT9", 
               "KIT", "EP300", "NF1", "COL1A1")
VlnPlot(combined_seurat, features = Marker.lt, ncol = 4)
FeaturePlot(combined_seurat, features = Marker.lt, ncol = 4)

DimPlot(combined_seurat, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "TissueType") + NoLegend()
DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "TissueType")
DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "SampleID" ) 
DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "seurat_clusters" ) 


# Load Seurat, ggplot2, and patchwork
library(Seurat)
library(patchwork)
library(ggplot2)

# 定義所有標記基因
all_marker_genes <- c(
  "KRT5", "KRT14", "TP63",       # KC_basal
  "KRT1", "KRT10", "IVL", "FLG", # KC_differentiating
  "PMEL", "DCT", "MITF",         # Melanocyte
  "COL6A3", "COL7A1", "MFAP5", "DPT", # Fibroblast
  "PDGFRB", "ACTA2", "MYH11",    # Pericyte/SMC
  "PECAM1", "VWF", "CLDN5",      # VEC
  "PROX1", "LYVE1", "PDPN",      # LEC
  "NCAM1", "KLRF1", "CD8A", "GZMB", "PRF1", # NK/Tcyt
  "CD3D", "CD3E", "CD69",        # Tcell_res
  "CD86", "IL1B", "NOS2",        # M1_like
  "CD163", "MRC1", "IL10",       # M2_like
  "CD14", "FCGR3A", "CD68",      # moMac
  "IL3RA", "GZMB", "CD1c", "CD11b", # DC_p
  "TPSAB1", "CPA3", "KIT"        # Mastcell
)

FeaturePlot(combined_seurat, features = all_marker_genes,  reduction = "umap", ncol = 6) 




# find markers for every cluster compared to all remaining cells, report only the positive
# ones
combined_seurat.markers <- FindAllMarkers(combined_seurat, only.pos = TRUE)
combined_seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


combined_seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10

DefaultAssay(combined_seurat) <- "RNA"
DoHeatmap(combined_seurat, features = top10$gene, assay = "integrated") + NoLegend()

library(Seurat)

# 繪製熱圖，使用自訂顏色範圍
DoHeatmap(combined_seurat, features = top10$gene, assay = "integrated") + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  NoLegend()


DoHeatmap(combined_seurat, features = top10$gene) + NoLegend()



#############################################################

# 假設您的Seurat物件名稱為 scRNA_data
# 定義每個Cluster的細胞類型名稱
celltype_annotations <- c(
  "Fibroblasts",        # Cluster 0
  "Differentiating Keratinocytes", # Cluster 1
  "Basal Keratinocytes", # Cluster 2
  "Dendritic Cells / Macrophages", # Cluster 3
  "Endothelial Cells",   # Cluster 4
  "Smooth Muscle Cells / Pericytes", # Cluster 5
  "Monocytes / Macrophages", # Cluster 6
  "T Cells",             # Cluster 7
  "NK Cells",            # Cluster 8
  "Mast Cells",          # Cluster 9
  "Melanocytes",         # Cluster 10
  "Lymphatic Endothelial Cells", # Cluster 11
  "Dendritic Cells",     # Cluster 12
  "Neutrophils",         # Cluster 13
  "Antigen-presenting Cells", # Cluster 14
  "Activated Fibroblasts", # Cluster 15
  "Pericytes",           # Cluster 16
  "Growth Factor-producing Fibroblasts", # Cluster 17
  "Smooth Muscle Cells"  # Cluster 18
)

# 將 seurat_clusters 轉換為數值型別，然後應用 celltype_annotations
combined_seurat$celltype <- factor(
  celltype_annotations[as.numeric(as.character(combined_seurat$seurat_clusters)) + 1], 
  levels = unique(celltype_annotations)
)




DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "celltype" )
# 檢查結果
head(scRNA_data[[]])

# 確保 ggplot2 套件已安裝並載入
if(!require('ggplot2')) {install.packages('ggplot2'); library(ggplot2)}

# 畫出直方圖
ggplot(combined_seurat@meta.data, aes(x = celltype, fill = TissueType)) +
  geom_bar(color = "black", size = 1.2) +  # 黑色粗框
  labs(x = "Cell type", y = "Count", fill = "Tissue Type") +
  theme_minimal() +
  scale_x_discrete(drop = FALSE) +  # 確保所有 cluster 類別顯示在 x 軸
  theme(
    axis.text = element_text(size = 14, face = "bold"),       # x, y 軸標籤加大加粗
    axis.title = element_text(size = 16, face = "bold"),      # x, y 軸標題加大加粗
    legend.title = element_text(size = 14, face = "bold"),    # 圖例標題加大加粗
    legend.text = element_text(size = 12, face = "bold")      # 圖例文字加大加粗
  )

