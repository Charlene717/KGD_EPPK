## Ref: https://erilu.github.io/bulk-rnaseq-analysis/#Differential_expression_analysis
## Ref: https://bioinformatics-core-shared-training.github.io/RNAseq-R/rna-seq-de.nb.html
## Ref: https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


##### Presetting #####
rm(list = ls()) # 清空環境

# 請先確保以下套件均已安裝
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if(!require("ggsignif")) install.packages("ggsignif"); library(ggsignif)
if(!require("readxl")) install.packages("readxl"); library(readxl)
if(!require("tibble")) install.packages("tibble"); library(tibble)

if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db"); library(org.Hs.eg.db)
if(!require("DESeq2")) BiocManager::install("DESeq2"); library(DESeq2)
if(!require("edgeR")) BiocManager::install("edgeR"); library(edgeR)

##### Set Parameters #####
# 改為多基因清單
Set_GeneName <- c("EPGN", "EGFR")
Set_DataType <- "GSE192832_Ichthyosis"

# 產生輸出資料夾(可自行更動)
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 10)
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))
Set_note <- paste0(Name_FileID, "_", Set_DataType)
Name_Export <- paste0(Name_FileID)
Name_ExportFolder <- paste0("Export_", Set_note)

if (!dir.exists(Name_ExportFolder)) dir.create(Name_ExportFolder)

##### Load Data #####
# 檔案路徑請依照實際情況修正
count_file <- "C:/Users/q2330/Dropbox/KGD_Lab/20241129_(PPK) Bulk RNA分析/GSE192832/GSE192832_raw_counts_GRCh38.p13_NCBI.tsv.gz"
annotation_xlsx <- "C:/Users/q2330/Dropbox/KGD_Lab/20241129_(PPK) Bulk RNA分析/GSE192832/GSE192832_Annotation.xlsx"

annotation_data <- readxl::read_excel(annotation_xlsx)
count_data <- read.table(gzfile(count_file), header = TRUE, row.names = 1, sep = "\t")

# 提取 row.names 作為 NCBI Gene IDs
gene_ids <- row.names(count_data)

# 基因名稱轉換 (ENTREZID -> SYMBOL)
gene_mapping <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys = gene_ids,
  keytype = "ENTREZID",
  columns = c("SYMBOL", "GENENAME")
)

gene_mapping <- gene_mapping %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

count_data <- as.data.frame(count_data) %>%
  rownames_to_column(var = "ENTREZID") %>%
  left_join(gene_mapping, by = "ENTREZID") %>%
  filter(!is.na(SYMBOL)) %>%
  column_to_rownames(var = "SYMBOL")

# 去除多餘欄位 (若原始檔最後多了幾欄非表達量的資訊)
# 這裡的 -1, -(ncol(count_data)) 若和你真實的資料欄位數對不上，請自行調整或移除。
# 如果 data 剛好已經正確，即無需以下範例操作，可註解掉。
count_data <- count_data[, c(-1, -(ncol(count_data)))]

##### Annotation 處理 #####
# 提取 subtype (NS, CIE, EI, LI, Control)
annotation_data <- annotation_data %>%
  mutate(subtype = gsub("subtype: ", "", `!Sample_characteristics_ch1...4`))

# 設定我們要關注的 subtype
subtypes <- c("NS", "CIE", "EI", "LI", "Control")

# 根據每個 subtype 找到對應的 SampleID
samples_list <- lapply(subtypes, function(stp) {
  annotation_data %>% 
    filter(subtype == !!stp) %>% 
    pull(GEONo)  # 或 pull(ID)，視實際欄位命名而定
})
names(samples_list) <- subtypes


##### Normalization #####
Set_NormM <- "DESeq2"  # 與 GSE54456/Urticaria 一致，也可改成 CPM / LogNor / None

if(Set_NormM == "DESeq2"){
  # 由於 DESeq2 需要 colData，其中必須包含設計公式中的欄位 (這裡是 ~ subtype)
  # 故先將 annotation_data 塞入
  # 注意 SampleID 與 count_data colnames 的對應，若不一致要想辦法對應
  # e.g. colnames(count_data) = GEONo ?
  
  # 若 count_data colnames 就是 GEONo，則：
  # colData: rownames 必須是樣本 ID (與 count_data 欄位一一對應)
  df_coldata <- annotation_data
  rownames(df_coldata) <- df_coldata$GEONo  # 假設 GEONo 與 colnames(count_data) 相符
  
  # 建立 DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = count_data, 
                                colData = df_coldata, 
                                design = ~ subtype)
  dds <- estimateSizeFactors(dds)
  # dds <-DESeq(dds)
  normalized_data <- counts(dds, normalized=TRUE)
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

##### 將表達資料分割成多基因清單 #####
gene_data_list <- list()

for (gene in Set_GeneName) {
  
  if (!gene %in% rownames(count_data)) {
    message(paste0("[Warning] ", gene, " not found in expression data! Skip."))
    next
  }
  
  # 取該基因的表達值
  gene_expr <- as.numeric(count_data[gene, ])
  # 建立一個小 data.frame: (SampleID, Expression)
  df_tmp <- data.frame(Expression = gene_expr, SampleID = colnames(count_data))
  
  # 接著根據 samples_list，將每個 subtype 內的樣本找出對應表達值
  # 或者直接用 annotation_data 去 merge 也可以
  # 這裡演示一下簡易做法：
  df_merged <- merge(df_tmp, annotation_data, by.x = "SampleID", by.y = "GEONo")
  
  # 這樣我們就有了 Expression 與 subtype
  df_merged$subtype <- factor(df_merged$subtype, levels = subtypes)
  
  # 存到 gene_data_list
  gene_data_list[[gene]] <- df_merged[, c("SampleID", "Expression", "subtype")]
}

##### 依照前例，多基因分別畫 BoxPlot #####
plot_list <- list()

# 設定比較組 (Control vs 其他四組)
comparisons <- list(
  c("Control", "LI"),
  c("Control", "EI"),
  c("Control", "CIE"),
  c("Control", "NS")
)

for (gene in names(gene_data_list)) {
  df_gene <- gene_data_list[[gene]]
  
  # 計算每個 subtype 有多少樣本
  sample_counts <- table(df_gene$subtype)
  
  Plot_Box <- ggplot(df_gene, aes(x = subtype, y = Expression, fill = subtype)) +
    geom_boxplot(outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.2, alpha = 0.7) +
    # 加上顯著性標記
    geom_signif(
      comparisons = comparisons,
      map_signif_level = TRUE,
      test = "wilcox.test",
      step_increase = 0.08
    ) +
    labs(
      title = paste(gene, "Expression in GSE192832 subtypes vs Control"),
      x = "Subtype",
      y = paste0(gene, " Expression"),
      caption = paste0(
        "NS: Netherton syndrome (n = ", sample_counts["NS"], "), ",
        "CIE: Congenital ichthyosiform erythroderma (n = ", sample_counts["CIE"], "), \n",
        "EI: Epidermolytic ichthyosis (n = ", sample_counts["EI"], "), ",
        "LI: Lamellar ichthyosis (n = ", sample_counts["LI"], "), ",
        "Control (n = ", sample_counts["Control"], ")"
      )

    ) +
    theme_minimal() +
    # 依喜好重新配置顏色
    scale_fill_manual(values = c(
      "NS"="#e85a5a","CIE"="#f78383","EI"="#fcb3b3",
      "LI"="#f5d3d3","Control"="#9999FF"
    )) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, face = "bold"),
      legend.position = "none",
      aspect.ratio = 1,
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  
  plot_list[[gene]] <- Plot_Box
}

# 輸出 BoxPlot 到 PDF
pdf(file = paste0(Name_ExportFolder, "/", Name_Export, "_BoxPlot.pdf"), width = 7, height = 5)
for (gene in names(plot_list)) {
  print(plot_list[[gene]])
}
dev.off()

##### 統計分析 (計算 log2FC、p-value...) #####
# 依照 GSE54456_Urticaria 的做法，定義函數
calculate_stats <- function(data, gene, group1, group2, test_method = "wilcox.test") {
  
  df2 <- data %>% filter(subtype %in% c(group1, group2))
  
  # 算平均值
  group_means <- df2 %>%
    group_by(subtype) %>%
    summarise(mean_expression = mean(Expression), .groups="drop") %>%
    arrange(match(subtype, c(group1, group2)))
  
  # FC = mean(group1) / mean(group2)
  # 此處以 "Control" vs "LI" 為例，要看你希望「 FC = group1 - group2 或 group2 - group1 」
  # 可自行調整順序
  FC <- group_means$mean_expression[2] / group_means$mean_expression[1]
  log2FC <- log2(FC)
  
  # 算 p-value
  x <- df2 %>% filter(subtype == group1) %>% pull(Expression)
  y <- df2 %>% filter(subtype == group2) %>% pull(Expression)
  
  if(test_method == "wilcox.test"){
    res <- wilcox.test(x, y)
  } else {
    res <- t.test(x, y)
  }
  
  data.frame(
    Gene = gene,
    Group1 = group1,
    Group2 = group2,
    log2FC = log2FC,
    FC = FC,
    p_value = res$p.value,
    statistic = res$statistic,
    test_method = test_method
  )
}

# 處理所有基因、所有比較
stats_results <- list()
for (gene in names(gene_data_list)) {
  df_gene <- gene_data_list[[gene]]
  
  for (comp in comparisons) {
    g1 <- comp[1]  # Control
    g2 <- comp[2]  # LI or EI or ...
    tmp <- calculate_stats(df_gene, gene, g1, g2, test_method = "wilcox.test")
    stats_results[[length(stats_results)+1]] <- tmp
  }
}

stats_results <- do.call(rbind, stats_results)
stats_results$p_adjusted <- p.adjust(stats_results$p_value, method = "BH")

# 可以額外加上資料集信息、疾病名稱等欄位
stats_results$Dataset <- "GSE192832"
stats_results$Disease <- "Ichthyosis" 
stats_results$DataType <- "RNA-seq"

# 儲存統計結果
write.csv(stats_results, paste0(Name_ExportFolder, "/", Name_Export, "_stats.csv"), row.names = FALSE)

##### 分布圖 (Distribution Plot) #####
# 將所有基因的資料合併為 long format
plot_data_long <- do.call(rbind, lapply(names(gene_data_list), function(gene) {
  tmp <- gene_data_list[[gene]]
  tmp$Gene <- gene
  tmp
}))

# 畫出 density 分布，並依 Gene 分面
Plot_Distrib <- ggplot(plot_data_long, aes(x = Expression, fill = subtype)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Gene, scales = "free", ncol = 2) +
  labs(
    title = "Expression Distribution by Marker",
    x = "Expression",
    y = "Density",
    fill = "Subtype"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c(
    "NS"="#e85a5a","CIE"="#f78383","EI"="#fcb3b3",
    "LI"="#f5d3d3","Control"="#9999FF"
  )) +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

pdf(file = paste0(Name_ExportFolder, "/", Name_Export, "_DistribPlot.pdf"), width = 8, height = 6)
print(Plot_Distrib)
dev.off()

##### (Optional) Export RData #####
save.image(file = paste0(Name_ExportFolder, "/", Name_Export, ".RData"))
