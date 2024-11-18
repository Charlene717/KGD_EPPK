##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
if(!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if(!require('Seurat')) install.packages('Seurat'); library(Seurat)
# remotes::install_version("Seurat", version = "4.3.0")



if(!require("scPred")) devtools::install_github("powellgenomicslab/scPred"); library(scPred)
# #open for Seurat Mult # trace("project_query", edit=TRUE) # layer = "data"
# trace("project_query", edit=TRUE) # GetAssayData(new, "data") -> # GetAssayData(new, "RNA")




load("D:/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/GSE202352_Normal_Batch effect correction_Temp.RData")
# load("C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/GSE202352_Normal_Batch effect correction_Temp.RData")

seuratObject_Sample <- combined_seurat

load("D:/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/GSE189889_Acral_Melanoma_Temp.RData")
# load("C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/GSE189889_Acral_Melanoma_Temp.RData")

seuratObject_Ref <- acral_melanoma
seuratObject_Ref_ori <- seuratObject_Ref

# # 移除 blueprint.main 為 NA 的細胞
# seuratObject_Ref_filtered <- seuratObject_Ref[, !is.na(seuratObject_Ref$blueprint.main)]
# 過濾掉 blueprint.main 為 "NA" 的細胞
seuratObject_Ref_filtered <- seuratObject_Ref[, seuratObject_Ref$blueprint.main != "NA"]

# 隨機抽取 5000 個細胞
set.seed(123)  # 設定隨機種子以確保結果可重現
selected_cells <- sample(colnames(seuratObject_Ref_filtered), 5000)

# 建立包含選定細胞的子集
seuratObject_Ref_filtered <- seuratObject_Ref_filtered[, selected_cells]


# 計算 blueprint.main 中各細胞類型的數量
cell_counts <- table(seuratObject_Ref_filtered$blueprint.main)

# 找出數量大於等於 5 的細胞類型
valid_types <- names(cell_counts[cell_counts >= 5])

# 過濾出 blueprint.main 屬於有效細胞類型的細胞
seuratObject_Ref <- seuratObject_Ref_filtered[, seuratObject_Ref_filtered$blueprint.main %in% valid_types]





seuratObject_Ref$blueprint.main %>% table()

seuratObject_Ref <- getFeatureSpace(seuratObject_Ref, "blueprint.main")   ## Get the feature space to train the classifiers
seuratObject_Ref <- trainModel(seuratObject_Ref) ## Train the model
# seuratObject_Ref <- trainModel(
#   object = seuratObject_Ref,          # Seurat 對象，需經過 getFeatureSpace
#   model = "svmRadial",                # 使用支援向量機（徑向基核）模型
#   preProcess = c("center", "scale"),  # 預處理選項：居中和標準化
#   resampleMethod = "cv",              # 使用交叉驗證
#   number = 3,                         # 減少交叉驗證摺數（例如設置為3）
#   seed = 66,                          # 固定隨機種子
#   tuneLength = 2,                     # 簡化參數網格，設置為2
#   metric = "Accuracy",                # 使用 Accuracy 作為評估指標來簡化訓練
#   returnData = FALSE,                 # 不返回訓練數據
#   savePredictions = "final",          # 保存最終預測
#   allowParallel = FALSE               # 不允許並行處理
# )

# DefaultAssay(seuratObject_Sample) <- "RNA"
# DefaultAssay(seuratObject_Ref) <- "RNA"
seuratObject_Sample[["data"]] <- seuratObject_Sample[["RNA"]]
seuratObject_Ref[["data"]] <- seuratObject_Ref[["RNA"]]
seuratObject_Sample <- scPredict(seuratObject_Sample, seuratObject_Ref) #, threshold = Set_scPredict_Thr)
DimPlot(seuratObject_Sample, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "scpred_prediction") + NoLegend()
DimPlot(seuratObject_Sample, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "scpred_prediction")
DimPlot(seuratObject_Sample, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "scpred_no_rejection") + NoLegend()
DimPlot(seuratObject_Sample, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "scpred_no_rejection")

colnames(seuratObject_Sample@meta.data)

DefaultAssay(seuratObject_Sample) <- "RNA"
Marker.lt <- c("EPGN", "EGFR","KRT1","KRT9",
               "KIT","EP300","NF1","COL1A1")
VlnPlot(seuratObject_Sample, features = Marker.lt,ncol = 4)

VlnPlot(seuratObject_Sample, features = Marker.lt[1:4],ncol = 2, group.by = "scpred_prediction")
VlnPlot(seuratObject_Sample, features = Marker.lt[5:8],ncol = 2, group.by = "scpred_prediction")


FeaturePlot(combined_seurat, features = Marker.lt, ncol = 4)

DimPlot(combined_seurat, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "GSM_ID") + NoLegend()
DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "GSM_ID")
DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "Tissue" ) 
DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "seurat_clusters" ) 
