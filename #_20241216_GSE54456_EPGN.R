##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)


# 檢查並載入必要的套件
# 確保必要的套件已安裝並載入
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("limma")) install.packages("limma"); library(limma)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if(!require("ggfortify")) install.packages("ggfortify"); library(ggfortify)

if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("EnhancedVolcano")) BiocManager::install("EnhancedVolcano"); library(EnhancedVolcano)
if(!require("ggsignif")) install.packages("ggsignif"); library(ggsignif)



#### Set Parameters ####
Target_Gene <- "EPGN"  # 可替換為其他基因名稱
Set_Dataset <- "GSE54456"

Set_NormM <- "DESeq2" # "DESeq2" # "None"
Set_PCA_3D <- TRUE

# Set export parameters
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 10) # Generate unique time-based ID
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))

Set_note <- paste0(Name_FileID,"_",Set_Dataset)


Name_Export <- paste0(Name_FileID)
Name_ExportFolder <- paste0("Export_",Set_note)

# Create export folder if it does not exist
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}


#### Load Data ####
# 讀取數據
count_file <- "C:/Users/q2330/Dropbox/KGD_Lab/#_數據分析_協助/Chitosan Microneedle/20241121 (鈺倫) 乾癬/Online_Dataset/GSE54456_raw_counts_GRCh38.p13_NCBI.tsv.gz"
annotation_file <- "C:/Users/q2330/Dropbox/KGD_Lab/#_數據分析_協助/Chitosan Microneedle/20241121 (鈺倫) 乾癬/Online_Dataset/GSE54456_Annotation_M.csv"

# 讀取 count 數據
count_data <- read.table(gzfile(count_file), header = TRUE, row.names = 1, sep = "\t")

# 讀取註解文件
annotation_data <- read.csv(annotation_file, header = TRUE)

# 設定分組因子
annotation_data <- annotation_data %>%
  mutate(Group = ifelse(Group == "normal_skin", "Normal", "Psoriasis"))

# 確保樣本名稱一致
annotation_data <- annotation_data[match(colnames(count_data), annotation_data$SampleID), ]


# 設置分組因子順序，確保 "Normal" 為參考組
group <- factor(annotation_data$Group, levels = c("Normal", "Psoriasis"))



#### PCA Analysis ####
# 過濾表達值全為零或變異性極低的基因
log2_data <- log2(count_data + 1)

# 計算基因的標準差
gene_sd <- apply(log2_data, 1, sd)

# 過濾標準差為零的基因
filtered_data <- log2_data[gene_sd > 0, ]

# PCA 分析
pca <- prcomp(t(filtered_data), scale. = TRUE)

# # 繪製 PCA 圖
# autoplot(pca, data = annotation_data, colour = 'Group', main = "PCA Plot")

# 提取 PCA 結果
pca_data <- as.data.frame(pca$x[, 1:2])  # 提取前兩個主成分
pca_data$ID <- rownames(pca_data)

# 合併 PCA 結果與 Annotation
merged_pca_data <- merge(pca_data, annotation_data, by.x = "ID", by.y = "SampleID")

# 繪製修訂版的 PCA 圖
library(ggplot2)
library(ggrepel)
# library(ggforce)

Plot_PCA <- ggplot(merged_pca_data, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 2, stroke = 1.5) +  # 設定點的大小，顏色顯示在圖例中
  geom_text_repel(aes(label = ID), size = 4, box.padding = 0.5, point.padding = 0.5) +  # 加大標註文字且不顯示在圖例中
  stat_ellipse(aes(color = Group), type = "norm", linetype = "dashed", size = 1) +  # 添加橢圓來圈出每個 cluster，不顯示圖例
  # geom_mark_hull(aes(fill = Group), concavity = 1, expand = 0.01, alpha = 0.2) +
  scale_color_manual(values = c("Psoriasis" = "#AE0000", "Normal" = "#6F00D2")) +  # 設定顏色
  theme_minimal(base_size = 6) +  # 設定整體字體大小
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),  # 添加黑色邊框
        axis.text = element_text(size = 16),  # 設定軸文字大小
        axis.title = element_text(size = 20), # 設定軸標題文字大小
        legend.text = element_text(size = 16), # 設定圖例文字大小
        legend.title = element_text(size = 20, face = "bold"),  # 設定圖例標題文字大小並加粗
        aspect.ratio = 1, # 固定 XY 軸長寬比例為 1:1
        legend.position = "right") +  # 圖例位置
  labs(x = paste0("PC 1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC 2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"))  # 設置圖例標題，只使用顏色來區分

# 顯示 PCA 圖
print(Plot_PCA)

pdf(paste0(Name_ExportFolder,"/", Name_Export,"_Plot_PCA.pdf"), width = 8, height = 10)

print(Plot_PCA)

dev.off()


####################################################################################
# 標準化資料 (Normalization)
if(Set_NormM == "DESeq2"){
  ## Ref: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = annotation_data, design = ~ Group)
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds)
  normalized_data <- counts(dds, normalized=TRUE)
  count_data <- counts(dds, normalized=TRUE)
  
}else if(Set_NormM == "LogNor"){
  # Log-normalization
  normalized_data <- log1p(count_data)
  
}else if(Set_NormM == "CPM"){
  # Counts per Million
  # normalized_data <- t(t(count_data) / colSums(count_data)) * 1e6
  normalized_data <- edgeR::cpm(count_data)
  count_data <- edgeR::cpm(count_data)
  # }else if(Set_NormM == "SCTransform"){
  #   if (!require("Seurat")) install.packages("Seurat"); library(Seurat)
  #   sct_data <- SCTransform(count_data, assay = "RNA")
  #   normalized_data <- GetAssayData(sct_data, slot = "data")
  
}else if(Set_NormM == "None"){
  
  normalized_data <- count_data
}





####################################################################################
# 設定分組因子
annotation_data <- annotation_data %>%
  mutate(Group = case_when(
    Group == "normal_skin" ~ "Normal",
    Group == "psoriasis" ~ "Psoriasis",
    TRUE ~ Group  # 保留原始值
  ))

group <- factor(annotation_data$Group, levels = c("Psoriasis", "Normal"))



if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {BiocManager::install("org.Hs.eg.db")}; library(org.Hs.eg.db) # 人類基因註釋資料庫

library(dplyr)        # 數據操作
library(tibble)       # 資料框轉換

# 提取 row.names 作為 NCBI Gene IDs
gene_ids <- row.names(count_data)

# 查詢基因名稱 (SYMBOL) 與描述 (GENENAME)
gene_mapping <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys = gene_ids,                 # 提供的 NCBI Gene IDs
  keytype = "ENTREZID",            # 指定輸入的基因 ID 類型
  columns = c("SYMBOL", "GENENAME")# 返回基因名稱和描述
)

# 檢查並移除重複的 SYMBOL，保留唯一值
gene_mapping <- gene_mapping %>%
  filter(!is.na(SYMBOL)) %>%         # 移除 SYMBOL 為 NA 的行
  distinct(SYMBOL, .keep_all = TRUE) # 保留唯一的 SYMBOL

# 將基因名稱合併到 count_data
# 確保 count_data 是 data.frame 格式
count_data <- as.data.frame(count_data)

# 合併 gene_mapping 到 count_data
count_data <- count_data %>%
  rownames_to_column(var = "ENTREZID") %>%   # 將 row.names 轉為一列
  left_join(gene_mapping, by = "ENTREZID") %>% # 左聯結基因名稱
  filter(!is.na(SYMBOL)) %>%                  # 移除無法匹配到基因名稱的行
  column_to_rownames(var = "SYMBOL")          # 使用 SYMBOL 作為新的 row.names



target_gene_expression <- count_data[Target_Gene, ]

# 整理 Normal 和 Psoriasis 組的目標基因表達值
normal_expression <- target_gene_expression[annotation_data$SampleID[annotation_data$Group == "Normal"]]
psoriasis_expression <- target_gene_expression[annotation_data$SampleID[annotation_data$Group == "Psoriasis"]]

# 統計分析 (t-test)
t_test_result <- t.test(as.numeric(normal_expression), as.numeric(psoriasis_expression), alternative = "two.sided")

# 可視化目標基因表達數據
data_to_plot <- data.frame(
  Expression = c(as.numeric(normal_expression), as.numeric(psoriasis_expression)),
  Group = factor(c(rep("Normal", length(normal_expression)), rep("Psoriasis", length(psoriasis_expression))),
                 levels = c("Psoriasis", "Normal"))
)

ggplot(data_to_plot, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_signif(
    comparisons = list(c("Psoriasis", "Normal")),
    map_signif_level = TRUE,
    y_position = max(data_to_plot$Expression) + 5,
    color = "black",
    textsize = 6,
    face = "bold",
    test = "wilcox.test"
  ) +
  labs(
    title = paste0(Target_Gene, " Expression in Psoriasis vs Normal (", Set_Dataset, ")"),
    x = "Group",
    y = paste0(Target_Gene, " Expression"),
    caption = paste0("Psoriasis: n = ", length(psoriasis_expression), "; Normal: n = ", length(normal_expression), ".")
  ) +
  scale_fill_manual(values = c("Psoriasis" = "#FF9999", "Normal" = "#9999FF")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    plot.caption = element_text(hjust = 0, size = 10),
    legend.position = "none",
    aspect.ratio = 1,
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
