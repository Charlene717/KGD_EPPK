# 檢查並載入必要套件
if (!require('DESeq2')) { install.packages('DESeq2'); library(DESeq2) }

# 提取數據
Marker.lt <- c("EPGN", "EGFR")
plot_data <- FetchData(seurat_object, vars = c(Marker.lt, "group"))

# 定義比較組別
group_comparisons <- list(
  c("AD_NON LESIONAL", "AD_LESIONAL"),
  c("LP_NON LESIONAL", "LP_LESIONAL"),
  c("Pso_NON LESIONAL", "Pso_LESIONAL")
)

# 構建 DESeq2 對象
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(seurat_object@assays$RNA@counts[Marker.lt, ]),
  colData = seurat_object@meta.data,
  design = ~ group
)

# 差異表達分析
dds <- DESeq(dds)  # 執行差異表達分析

# 初始化結果儲存
all_results <- data.frame()

# 計算每組比較的 logFC 和 p 值
for (comparison in group_comparisons) {
  contrast_name <- paste(comparison, collapse = " vs ")
  
  # 獲取差異分析結果
  res <- results(dds, contrast = c("group", comparison[2], comparison[1]))
  
  # 選取感興趣的基因
  res_filtered <- as.data.frame(res[rownames(res) %in% Marker.lt, ])
  
  # 添加比較組別標籤
  res_filtered$Group1 <- comparison[1]
  res_filtered$Group2 <- comparison[2]
  
  # 添加比較名稱
  res_filtered$Comparison <- contrast_name
  
  # 合併結果
  all_results <- rbind(all_results, res_filtered)
}

# 校正 p 值
all_results$p_adjusted <- p.adjust(all_results$pvalue, method = "BH")

# 查看結果
print(all_results)

# 將結果保存到文件（可選）
write.csv(all_results, "deseq2_group_comparisons_results.csv", row.names = TRUE)
