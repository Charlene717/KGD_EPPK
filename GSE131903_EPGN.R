# 設定目標基因名稱
Target_Gene <- "EPGN"  # 可自由替換為其他基因名稱

# 確保必要的套件是否已安裝
if(!require("readxl")) {install.packages("readxl"); library(readxl)}
if(!require("dplyr")) {install.packages("dplyr"); library(dplyr)}
if(!require("biomaRt")) {install.packages("biomaRt"); library(biomaRt)}
if(!require("ggplot2")) {install.packages("ggplot2"); library(ggplot2)}
if(!require("ggpubr")) {install.packages("ggpubr"); library(ggpubr)}

# 設定檔案路徑
file_path <- "C:/Users/q2330/Dropbox/KGD_Lab/20241129_(PPK) Bulk RNA分析/GSE131903/GSE131903_HI_norm_counts.xlsx"

# 讀取資料
data <- read_excel(file_path)

# 調整列名
colnames(data)[1] <- "Ensembl_ID"

# 基因編碼轉換為標準名稱（使用 biomaRt）
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "uswest")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = data$Ensembl_ID,
  mart = ensembl
)

# 附加標準基因名稱到數據框
data <- data %>% 
  left_join(gene_map, by = c("Ensembl_ID" = "ensembl_gene_id"))

# 篩選目標基因數據
target_gene_data <- data %>% filter(external_gene_name == Target_Gene)

if (nrow(target_gene_data) == 0) {
  stop(paste("Gene", Target_Gene, "not found in the dataset. Please check the gene name."))
}

# 確定組別（從欄位名稱）
ns_samples <- grep("^NS", colnames(target_gene_data), value = TRUE)
hi_samples <- grep("^HI", colnames(target_gene_data), value = TRUE)

# 擷取數據並轉換為數值（對數轉換）
target_gene_ns <- log2(as.numeric(as.matrix(target_gene_data[, ns_samples])) + 1)
target_gene_hi <- log2(as.numeric(as.matrix(target_gene_data[, hi_samples])) + 1)

# t 檢定
t_test_result <- t.test(target_gene_hi, target_gene_ns, alternative = "greater")

# 數據準備（結合兩組數據）
target_gene_combined <- data.frame(
  Expression = c(target_gene_hi, target_gene_ns),
  Group = c(rep("HI", length(target_gene_hi)), rep("NS", length(target_gene_ns)))
)

# 新增組別全名和樣本數
hi_sample_count <- length(target_gene_hi)
ns_sample_count <- length(target_gene_ns)
legend_text <- paste0(
  "Harlequin Ichthyosis (HI, n = ", hi_sample_count, 
  ") and Normal Skin (NS, n = ", ns_sample_count, 
  ")."
)

# 計算比較框線位置
y_max <- max(target_gene_combined$Expression) + 1

# 繪製圖表
ggplot(target_gene_combined, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_signif(
    comparisons = list(c("HI", "NS")),
    map_signif_level = TRUE,
    y_position = y_max,
    color = "black",
    textsize = 6,
    face = "bold",
    test = "t.test"
  ) +
  labs(
    title = paste(Target_Gene, "Expression in HI vs NS (GSE131903)"),
    x = "Group",
    y = "Log2 Transformed Expression",
    caption = legend_text
  ) +
  scale_fill_manual(values = c("HI" = "#FF9999", "NS" = "#9999FF")) +
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
