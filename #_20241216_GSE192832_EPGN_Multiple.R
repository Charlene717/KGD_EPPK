# 設定目標基因名稱
Target_Gene <- "EPGN"  # 可自由替換成其他基因名稱

# 確保必要套件已安裝
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if(!require("ggsignif")) install.packages("ggsignif"); library(ggsignif)
if(!require("readxl")) install.packages("readxl"); library(readxl)

# 讀取表達數據和註解數據
count_file <- "C:/Users/q2330/Dropbox/KGD_Lab/20241129_(PPK) Bulk RNA分析/GSE192832/GSE192832_raw_counts_GRCh38.p13_NCBI.tsv.gz"
annotation_data <- readxl::read_excel("C:/Users/q2330/Dropbox/KGD_Lab/20241129_(PPK) Bulk RNA分析/GSE192832/GSE192832_Annotation.xlsx")


# 讀取 count 數據
count_data <- read.table(gzfile(count_file), header = TRUE, row.names = 1, sep = "\t")


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

count_data <- count_data[,c(-1,-ncol(count_data))]

####################################################################################
Set_NormM <- "DESeq2"
# 標準化資料 (Normalization)
if(Set_NormM == "DESeq2"){
  ## Ref: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = annotation_data, design = ~ subtype)
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





# 提取樣本的 subtype 信息
annotation_data <- annotation_data %>%
  mutate(subtype = gsub("subtype: ", "", `!Sample_characteristics_ch1...4`))

# 篩選出需要的樣本 ID
subtypes <- c("NS", "CIE", "EI", "LI", "Control")
samples_list <- lapply(subtypes, function(subtype) {
  annotation_data %>%
    filter(subtype == !!subtype) %>%
    pull(GEONo)
    # pull(ID) %>%
    # paste0("X", .)  # 調整為表達數據的列名格式
})
names(samples_list) <- subtypes

# 檢查目標基因是否在表達數據中
if (Target_Gene %in% rownames(count_data)) {
  
  # 提取目標基因表達數據
  gene_expression <- count_data[Target_Gene, ]
  
  # 整理各 subtype 和 Control 組的目標基因表達值
  data_to_plot <- do.call(rbind, lapply(names(samples_list), function(subtype) {
    data.frame(
      Expression = as.numeric(gene_expression[samples_list[[subtype]]]),
      Group = subtype
    )
  }))
  
  # 計算樣本數
  sample_counts <- sapply(samples_list, length)
  
  # 進行統計分析 (t-test) 並標記顯著性
  comparisons <- list(
    c("Control", "LI"),
    c("Control", "EI"),
    c("Control", "CIE"),
    c("Control", "NS")
  )
  
  # 繪製圖表
  ggplot(data_to_plot, aes(x = factor(Group, levels = c("NS", "CIE", "EI", "LI", "Control")), y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.2, alpha = 0.7) +
    geom_signif(
      comparisons = comparisons,
      map_signif_level = TRUE,
      test = "wilcox.test",
      step_increase = 0.1
    ) +
    labs(
      title = paste(Target_Gene, "Expression in LS Subtypes vs Control (GSE192832)"),
      x = "Group",
      y = paste(Target_Gene, "Expression"),
      caption = paste0(
        "NS: Netherton syndrome (n = ", sample_counts["NS"], "), ",
        "CIE: congenital ichthyosiform erythroderma (n = ", sample_counts["CIE"], "), \n",
        "EI: epidermolytic ichthyosis (n = ", sample_counts["EI"], "), ",
        "LI: lamellar ichthyosis (n = ", sample_counts["LI"], "), ",
        "Control (n = ", sample_counts["Control"], ")"
      )
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("LI" = "#f5d3d3", "EI" = "#fcb3b3", "CIE" = "#f78383", "NS" = "#e85a5a", "Control" = "#9999FF")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 14, face = "bold"),
      plot.caption = element_text(hjust = 0, size = 10),
      legend.position = "none",
      aspect.ratio = 1,
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
} else {
  message("Gene ", Target_Gene, " is not found in the expression data. Please check the gene name.")
}
