##### Presetting #####
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Package #####
if (!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if (!require("limma")) install.packages("limma"); library(limma)
if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if (!require("ggfortify")) install.packages("ggfortify"); library(ggfortify)
if (!require("ggsignif")) install.packages("ggsignif"); library(ggsignif)
if (!require("EnhancedVolcano")) BiocManager::install("EnhancedVolcano"); library(EnhancedVolcano)

# Bioconductor / Annotation
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db"); library(org.Hs.eg.db)

# 其他分析工具
if (!require("DESeq2")) BiocManager::install("DESeq2"); library(DESeq2)
if (!require("edgeR")) BiocManager::install("edgeR"); library(edgeR)

##### Set Parameters #####
Set_GeneName <- c("EPGN", "EGFR")  # 多個基因名稱
Set_Dataset <- "GSE54456_Psoriasis"

# 若需要可自行切換其他 Normalization 方法，如： "LogNor", "CPM", "None" 等
Set_NormM <- "DESeq2" 
Set_PCA_3D <- TRUE

# Set export parameters
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 10) # Generate unique time-based ID
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))
Set_note <- paste0(Name_FileID,"_", Set_Dataset)

Name_Export <- paste0(Name_FileID)
Name_ExportFolder <- paste0("Export_", Set_note)

# Create export folder if it does not exist
if (!dir.exists(Name_ExportFolder)) { dir.create(Name_ExportFolder) }

##### Load Data #####
# 請依實際檔案路徑做調整
count_file <- "C:/Users/q2330/Dropbox/KGD_Lab/#_數據分析_協助/Chitosan Microneedle/20241121 (鈺倫) 乾癬/Online_Dataset/GSE54456_raw_counts_GRCh38.p13_NCBI.tsv.gz"
annotation_file <- "C:/Users/q2330/Dropbox/KGD_Lab/#_數據分析_協助/Chitosan Microneedle/20241121 (鈺倫) 乾癬/Online_Dataset/GSE54456_Annotation_M.csv"

# 讀取 count 數據
count_data <- read.table(gzfile(count_file), header = TRUE, row.names = 1, sep = "\t")

# 讀取註解文件
annotation_data <- read.csv(annotation_file, header = TRUE)

# 確保樣本名稱一致 & 分組標記
annotation_data <- annotation_data[match(colnames(count_data), annotation_data$SampleID), ]
annotation_data <- annotation_data %>%
  mutate(Group = ifelse(Group == "normal_skin", "Normal", "Psoriasis"))

# 設定分組因子 (可視需求改變順序：若想要 Normal 當參考組，就把 Normal 放前面)
group <- factor(annotation_data$Group, levels = c("Normal", "Psoriasis"))

##### Quick PCA (Optional) #####
# 過濾表達值全為零或變異性極低的基因
log2_data <- log2(count_data + 1)
gene_sd <- apply(log2_data, 1, sd)
filtered_data <- log2_data[gene_sd > 0, ]

# PCA
pca <- prcomp(t(filtered_data), scale. = TRUE)

# 取 PCA 前兩軸
pca_data <- as.data.frame(pca$x[, 1:2])  
pca_data$ID <- rownames(pca_data)
merged_pca_data <- merge(pca_data, annotation_data, by.x = "ID", by.y = "SampleID")

Plot_PCA <- ggplot(merged_pca_data, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 2, stroke = 1.5) +
  ggrepel::geom_text_repel(aes(label = ID), size = 4, box.padding = 0.5, point.padding = 0.5) +
  stat_ellipse(aes(color = Group), type = "norm", linetype = "dashed", size = 1) +
  scale_color_manual(values = c("Psoriasis" = "#AE0000", "Normal" = "#6F00D2")) +
  theme_minimal(base_size = 6) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20, face = "bold"),
    aspect.ratio = 1,
    legend.position = "right"
  ) +
  labs(
    x = paste0("PC 1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC 2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")
  )

print(Plot_PCA)
pdf(paste0(Name_ExportFolder,"/", Name_Export,"_Plot_PCA.pdf"), width = 8, height = 10)
print(Plot_PCA)
dev.off()

##### Normalization #####
if(Set_NormM == "DESeq2"){
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = annotation_data, design = ~ Group)
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds)
  normalized_data <- counts(dds, normalized = TRUE)
  count_data <- normalized_data
  
} else if(Set_NormM == "LogNor"){
  normalized_data <- log1p(count_data)
  count_data <- normalized_data
  
} else if(Set_NormM == "CPM"){
  normalized_data <- edgeR::cpm(count_data)
  count_data <- normalized_data
  
} else if(Set_NormM == "None"){
  normalized_data <- count_data
}

##### 基因名稱轉換 (ENTREZID -> SYMBOL) #####
# 以 row.names 為 ENTREZID
gene_ids <- row.names(count_data)

gene_mapping <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys = gene_ids,
  keytype = "ENTREZID",
  columns = c("SYMBOL", "GENENAME")
)

# 移除重複與 NA
gene_mapping <- gene_mapping %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

# 整理 count_data
count_data <- as.data.frame(count_data) %>%
  tibble::rownames_to_column(var = "ENTREZID") %>%
  left_join(gene_mapping, by = "ENTREZID") %>%
  filter(!is.na(SYMBOL)) %>%
  tibble::column_to_rownames(var = "SYMBOL")

##### 將表達資料分割成多基因清單 (類似 GSE57178_Urticaria 的做法) #####
# 與之前 GSE57178 的寫法一致，逐基因取資料
gene_data_list <- list()

for (gene in Set_GeneName) {
  # 如果該基因不存在，跳過
  if (!gene %in% rownames(count_data)) {
    message(paste0(gene, " not found in expression data!"))
    next
  }
  
  # 取該基因的表達值
  gene_expr <- as.numeric(count_data[gene, ])
  
  # 與 annotation_data 合併
  df_tmp <- data.frame(
    Expression = gene_expr,
    SampleID = colnames(count_data)
  )
  
  # 匹配分組
  df_merged <- merge(df_tmp, annotation_data, by.x = "SampleID", by.y = "SampleID")
  # 確保 Group 像之前一樣 (Normal, Psoriasis)
  df_merged$Group <- factor(df_merged$Group, levels = c("Normal", "Psoriasis"))
  
  # 存到 list
  gene_data_list[[gene]] <- df_merged[, c("SampleID", "Expression", "Group")]
}


##### 依照 GSE57178_Urticaria 的方式，分別為每個基因畫 BoxPlot #####
plot_list <- list()

for (gene in names(gene_data_list)) {
  data <- gene_data_list[[gene]]
  
  # 設定因子順序 (左到右: Psoriasis -> Normal)
  data$Group <- factor(data$Group, levels = c("Psoriasis", "Normal"))
  
  # 比較組別同樣設定為 Psoriasis 在前、Normal 在後
  comparisons <- list(c("Psoriasis", "Normal"))
  
  Plot_Box <- ggplot(data, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    geom_signif(
      comparisons = comparisons,
      map_signif_level = TRUE,
      test = "wilcox.test",
      step_increase = 0.2
    ) +
    labs(
      title = paste0(gene, " Expression in Psoriasis vs Normal"),
      x = "Group",
      y = paste0("Expression (", gene, ")")
    ) +
    scale_fill_manual(values = c("Psoriasis" = "#FF9999", "Normal" = "#9999FF")) +
    # 在這裡透過 scale_x_discrete() 指定 x 軸顯示順序
    scale_x_discrete(limits = c("Psoriasis", "Normal")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      legend.position = "none",
      aspect.ratio = 1,
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  
  plot_list[[gene]] <- Plot_Box
  print(Plot_Box)
}


# 輸出所有基因的 BoxPlot
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_BoxPlot.pdf"), width = 6, height = 6)
for (gene in names(plot_list)) {
  print(plot_list[[gene]])
}
dev.off()

##### Statistical Analysis (log2FC, p-value, etc.) #####
# 定義比較組
group_comparisons <- list(c("Psoriasis", "Normal"))

calculate_stats <- function(data, marker, groups, test_method = "wilcox.test") {
  # 僅保留有興趣的兩組
  df2 <- data %>% dplyr::filter(Group %in% groups)
  
  # 計算兩組的均值
  group_means <- df2 %>%
    group_by(Group) %>%
    summarize(mean_expression = mean(Expression), .groups = 'drop') %>%
    arrange(match(Group, groups))
  
  # FC = mean(Psoriasis) / mean(Normal)，此處的前後順序可自行調整
  FC <- group_means$mean_expression[1] / group_means$mean_expression[2]
  log2FC <- log2(FC)
  
  # p-value
  group1_values <- df2 %>% filter(Group == groups[1]) %>% pull(Expression)
  group2_values <- df2 %>% filter(Group == groups[2]) %>% pull(Expression)
  
  if (test_method == "t.test") {
    test_result <- t.test(group1_values, group2_values)
  } else if (test_method == "wilcox.test") {
    test_result <- wilcox.test(group1_values, group2_values)
  } else {
    stop("Unsupported test method")
  }
  
  p_value <- test_result$p.value
  statistic <- test_result$statistic
  
  data.frame(
    Marker = marker,
    Group1 = groups[1],
    Group2 = groups[2],
    log2FC = log2FC,
    FC = FC,
    p_value = p_value,
    statistic = statistic,
    test = test_method
  )
}

stats_results <- lapply(names(gene_data_list), function(marker) {
  do.call(rbind, lapply(group_comparisons, function(comp) {
    calculate_stats(gene_data_list[[marker]], marker, comp, test_method = "wilcox.test")
  }))
}) %>%
  do.call(rbind, .)

# 多重比較校正
stats_results$p_adjusted <- p.adjust(stats_results$p_value, method = "BH")

# 加入與 GSE57178_Urticaria 相同的描述欄位
stats_results$Disease <- "Psoriasis"
stats_results$Dataset <- "GSE54456"
stats_results$DataType <- "RNA-seq"

# 重新排序 rownames
rownames(stats_results) <- seq_len(nrow(stats_results))
print(stats_results)

# 輸出統計結果
write.csv(stats_results, paste0(Name_ExportFolder, "/", Name_Export, "_stats.csv"), row.names = TRUE)

##### Distribution Plot (多基因) #####
# 將數據轉長格式
plot_data_long <- do.call(rbind, lapply(names(gene_data_list), function(marker) {
  tmp <- gene_data_list[[marker]]
  tmp$Marker <- marker
  tmp
}))

# 根據需求，可用 log2(Expression+1) 再畫，也可直接 Expression
Plot_Distrib <- ggplot(plot_data_long, aes(x = Expression, fill = Group)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ Marker, scales = "free", ncol = 2) +
  labs(
    title = "Expression Distribution by Marker",
    x = "Expression",
    y = "Density",
    fill = "Group"
  ) +
  theme_minimal(base_size = 15) +
  scale_fill_manual(values = c("Psoriasis" = "#FF9999", "Normal" = "#9999FF")) +
  theme(
    strip.text = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

pdf(paste0(Name_ExportFolder, "/", Name_Export, "_DistribPlot.pdf"), width = 10, height = 6)
print(Plot_Distrib)
dev.off()

##### Export RData for record #####
save.image(paste0(Name_ExportFolder, "/", Name_Export, ".RData"))
