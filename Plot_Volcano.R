##### 火山圖 (Volcano Plot) #####

# 確保已安裝並載入必要套件
if(!require("ggrepel")) install.packages("ggrepel"); library(ggrepel)

# 執行 DESeq2 完整分析流程
dds <- DESeq(dds)

# 獲取差異表達結果
res <- results(dds)

# 將結果轉換為 data.frame 並加入基因名稱
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)

# 準備火山圖資料
volcano_data <- res_df %>%
  mutate(
    log10_p = -log10(pvalue),
    Highlight = ifelse(Gene %in% Set_GeneName, "Yes", "No")
  ) %>%
  filter(!is.na(log2FoldChange) & !is.na(pvalue))

# 繪製火山圖，確保Set_GeneName基因位於上層
Plot_Volcano <- ggplot() +
  # 所有非關注基因
  geom_point(data = subset(volcano_data, Highlight == "No"), 
             aes(x = log2FoldChange, y = log10_p), 
             color = "grey", alpha = 0.6) +
  # 關注基因
  geom_point(data = subset(volcano_data, Highlight == "Yes"), 
             aes(x = log2FoldChange, y = log10_p), 
             color = "red", alpha = 0.8) +
  # 標示基因名稱
  geom_text_repel(
    data = subset(volcano_data, Highlight == "Yes"),
    aes(x = log2FoldChange, y = log10_p, label = Gene),
    size = 5,
    box.padding = 0.3,
    point.padding = 0.5,
    max.overlaps = Inf
  ) +
  # 閾值線
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  # 標籤與主題
  labs(
    title = "Volcano Plot of Differential Expression",
    x = "Log2 Fold Change",
    y = "-Log10(p-value)",
    color = "Gene Highlight"
  ) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

Plot_Volcano

# 儲存火山圖到 PDF
pdf(file = paste0(Name_ExportFolder, "/", Name_Export, "_VolcanoPlot.pdf"), width = 8, height = 6)
print(Plot_Volcano)
dev.off()




# 繪製火山圖
Plot_Volcano_O <- ggplot(volcano_data, aes(x = log2FoldChange, y = log10_p)) +
  geom_point(aes(color = Highlight), alpha = 0.6) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = subset(volcano_data, Highlight == "Yes"),
    aes(label = Gene),
    size = 5,
    box.padding = 0.3,
    point.padding = 0.5,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano Plot of Differential Expression",
    x = "Log2 Fold Change",
    y = "-Log10(p-value)",
    color = "Gene Highlight"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
Plot_Volcano_O
