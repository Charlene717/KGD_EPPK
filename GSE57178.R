# 確保所需的套件已安裝
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("GEOquery")) BiocManager::install("GEOquery")
if (!require("limma")) BiocManager::install("limma")

# 加載必要的套件
library(GEOquery)
library(limma)

# 加載 GEO 數據
gset <- getGEO("GSE57178", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# 提取表達矩陣
ex <- exprs(gset)

# 提取註釋數據
platform <- getGEO("GPL6244", AnnotGPL = TRUE)
annot <- Table(platform)

# 確認基因名稱列是否正確 (例如 "Gene symbol")
epgn_row <- annot$ID[annot$`Gene symbol` == "EPGN"]

# 如果找不到，返回錯誤
if (length(epgn_row) == 0) stop("EPGN not found in annotation data!")

# 確認 ex 的行名和註釋的 ID 是否一致
rownames_ex <- rownames(ex)
matched_id <- intersect(epgn_row, rownames_ex)

if (length(matched_id) == 0) stop("No matching IDs found for EPGN in expression data!")

# 提取 EPGN 的表達數據
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
        main = "EPGN Expression",
        xlab = "Group",
        ylab = "Expression Level (log2)",
        col = c("pink","lightblue", "lightgreen"),
        border = "black", notch = FALSE)

# 添加統計顯著性檢驗 (t檢驗)
urticaria <- data$Expression[data$Group == "Urticaria"]
normal_skin <- data$Expression[data$Group == "Normal Skin"]
p_value <- t.test(urticaria, normal_skin)$p.value
legend("topright", legend = paste("p-value =", signif(p_value, 3)), bty = "n")
