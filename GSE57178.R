##### Presetting #####
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Package #####
if(!require('Seurat')) {install.packages('Seurat'); library(Seurat)}
if(!require('tidyverse')) {install.packages('tidyverse'); library(tidyverse)}
if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("GEOquery")) {BiocManager::install("GEOquery"); library(GEOquery)}
if (!require("limma")) {BiocManager::install("limma"); library(limma)}


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

Set_GeneName <- "EPGN"
# 確認基因名稱列是否正確 (例如 "Gene symbol")
epgn_row <- annot$ID[annot$`Gene symbol` == Set_GeneName]

# 如果找不到，返回錯誤
if (length(epgn_row) == 0) stop(paste0(Set_GeneName," not found in annotation data!"))

# 確認 ex 的行名和註釋的 ID 是否一致
rownames_ex <- rownames(ex)
matched_id <- intersect(epgn_row, rownames_ex)

if (length(matched_id) == 0) stop(paste0("No matching IDs found for ",Set_GeneName," in expression data!"))

# 提取 Set_GeneName 的表達數據
epgn_expr <- ex[matched_id, ]

# 如果有多個探針，取平均值
if (is.matrix(epgn_expr)) {
  epgn_expr <- rowMeans(epgn_expr)
}

# 獲取樣本分組信息
pdata <- pData(gset)
groups <- pdata$`characteristics_ch1`  # 假設分組信息在該列
groups <- gsub(".*condition: ", "", groups)  # 清理分組標籤

# 構建數據框
data <- data.frame(Expression = as.numeric(epgn_expr), Group = groups)

# 畫盒鬚圖
boxplot(Expression ~ Group, data = data,
        main = paste0(Set_GeneName, " Expression"),
        xlab = "Group",
        ylab = "Expression Level (log2)",
        col = c("pink","lightblue", "lightgreen"),
        border = "black", notch = FALSE)


######################

# 確保所需的套件已安裝
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggsignif")) install.packages("ggsignif")

# 加載所需套件
library(ggplot2)
library(ggsignif)

# 整理數據，縮寫組名
data$Group <- gsub("tissue: lesional skin", "Lesional", data$Group)
data$Group <- gsub("tissue: non-lesional skin", "Non-lesional", data$Group)
data$Group <- gsub("tissue: normal skin", "Normal", data$Group)

# 設定組別順序
data$Group <- factor(data$Group, levels = c("Lesional", "Non-lesional", "Normal"))

# 計算樣本數
sample_counts <- table(data$Group)

# 定義比較組
comparisons <- list(
  c("Lesional", "Non-lesional"),
  c("Lesional", "Normal")
)

# 繪製盒鬚圖
ggplot(data, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_signif(
    comparisons = comparisons,
    map_signif_level = TRUE,
    test = "t.test", # 可以改為 "t.test"、"wilcox.test" 根據需求
    step_increase = 0.25
  ) +
  labs(
    title = paste0(Set_GeneName," Expression in Skin Groups"),
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
  )+
  scale_y_continuous(limits = c(0, NA)) # 明確設置 y 軸從 0 開始

