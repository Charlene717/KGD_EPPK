##### Presetting #####
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Package #####
if(!require('Seurat')) {install.packages('Seurat'); library(Seurat)}
if(!require('tidyverse')) {install.packages('tidyverse'); library(tidyverse)}
if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}
if(!require('ggplot2')) {install.packages('ggplot2'); library(ggplot2)}
if(!require('ggsignif')) {install.packages('ggsignif'); library(ggsignif)}
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("GEOquery")) {BiocManager::install("GEOquery"); library(GEOquery)}
if (!require("limma")) {BiocManager::install("limma"); library(limma)}

##### Set Parameter #####
Set_GeneName <- c("EPGN", "EGFR") # 多個基因名稱
Set_DataType <- "GSE57178_Urticaria"

# Set export parameters
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 10) # Generate unique time-based ID
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))
Set_note <- paste0(Name_FileID, "_", Set_DataType)
Name_Export <- paste0(Name_FileID)
Name_ExportFolder <- paste0("Export_", Set_note)
# Create export folder if it does not exist
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}

##### Load Dataset #####
# 加載 GEO 數據
gset <- getGEO("GSE57178", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# 提取表達矩陣
ex <- exprs(gset)

# 提取註釋數據
platform <- getGEO("GPL6244", AnnotGPL = TRUE)
annot <- Table(platform)

##### Data Processing #####
# 初始化存儲表達數據的列表
gene_data_list <- list()

for (gene in Set_GeneName) {
  # 確認基因名稱列是否正確 (例如 "Gene symbol")
  gene_row <- annot$ID[annot$`Gene symbol` == gene]
  
  # 如果找不到，跳過此基因
  if (length(gene_row) == 0) {
    message(paste0(gene, " not found in annotation data!"))
    next
  }
  
  # 確認 ex 的行名和註釋的 ID 是否一致
  rownames_ex <- rownames(ex)
  matched_id <- intersect(gene_row, rownames_ex)
  
  if (length(matched_id) == 0) {
    message(paste0("No matching IDs found for ", gene, " in expression data!"))
    next
  }
  
  # 提取基因的表達數據
  gene_expr <- ex[matched_id, ]
  
  # 如果有多個探針，取平均值
  if (is.matrix(gene_expr)) {
    gene_expr <- rowMeans(gene_expr)
  }
  
  # 獲取樣本分組信息
  pdata <- pData(gset)
  groups <- pdata$`characteristics_ch1`  # 假設分組信息在該列
  groups <- gsub(".*condition: ", "", groups)  # 清理分組標籤
  
  # 構建數據框
  data <- data.frame(Expression = as.numeric(gene_expr), Group = groups)
  
  # 整理數據，縮寫組名
  data$Group <- gsub("tissue: lesional skin", "Lesional", data$Group)
  data$Group <- gsub("tissue: non-lesional skin", "Non-lesional", data$Group)
  data$Group <- gsub("tissue: normal skin", "Normal", data$Group)
  
  # 設定組別順序
  data$Group <- factor(data$Group, levels = c("Lesional", "Non-lesional", "Normal"))
  
  # 保存數據到列表中
  gene_data_list[[gene]] <- data
}

##### Visualization #####
# 初始化存儲圖形的列表
plot_list <- list()

# 開始生成圖形
for (gene in names(gene_data_list)) {
  data <- gene_data_list[[gene]]
  
  # 計算樣本數
  sample_counts <- table(data$Group)
  
  # 定義比較組
  comparisons <- list(
    c("Lesional", "Non-lesional"),
    c("Lesional", "Normal")
  )
  
  # 繪製盒鬚圖
  Plot_Box <- ggplot(data, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.2, alpha = 0.7) +
    geom_signif(
      comparisons = comparisons,
      map_signif_level = TRUE,
      test = "t.test", # 可以改為 "t.test" 或 "wilcox.test" 根據需求
      step_increase = 0.25
    ) +
    labs(
      title = paste0(gene, " Expression in Skin Groups"),
      x = "Group",
      y = "Expression Level (log2)",
      caption = paste0(
        "Lesional (n = ", sample_counts["Lesional"], "), ",
        "Non-lesional (n = ", sample_counts["Non-lesional"], "), ",
        "Normal (n = ", sample_counts["Normal"], ")"
      )
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("Normal" = "#9999FF", "Non-lesional" = "#f5d3d3", "Lesional" = "#e85a5a")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      plot.caption = element_text(hjust = 0, size = 10),
      legend.position = "none",
      aspect.ratio = 1,
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +
    scale_y_continuous(limits = c(0, NA)) # 明確設置 y 軸從 0 開始
  
  # 保存圖形到列表
  plot_list[[gene]] <- Plot_Box
}

# 將圖形輸出到 PDF
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_BoxPlot.pdf"), width = 6, height = 6)
for (gene in names(plot_list)) {
  print(plot_list[[gene]])
}
dev.off()

#### Export ####
## Export RData
save.image(paste0(Name_ExportFolder, "/", Name_Export, ".RData"))
