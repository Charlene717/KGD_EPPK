# 確保所需的套件已安裝
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("GEOquery")) BiocManager::install("GEOquery")
if (!require("limma")) BiocManager::install("limma")

# 加載必要的套件
library(GEOquery)
library(limma)

# 加載數據集
gset <- getGEO("GSE57178", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# 提取表達矩陣
ex <- exprs(gset)

# 檢查是否需要 log2 轉換
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  ex <- log2(ex)
}

# 獲取樣本的分組信息
pdata <- pData(gset)
groups <- pdata$`characteristics_ch1`  # 假設分組信息在該列
# 確定組名稱 (例如 "Urticaria" 和 "Normal Skin")
groups <- gsub(".*condition: ", "", groups)  # 根據數據清理分組標籤

# 獲取 EPGN 的表達數據
epgn_expr <- ex["EPGN", ]  # 假設行名包含基因名稱

# 合併分組與表達數據
data <- data.frame(Expression = epgn_expr, Group = groups)

# 畫盒鬚圖
boxplot(Expression ~ Group, data = data, 
        main = "EPGN Expression in Urticaria vs Normal Skin",
        xlab = "Group",
        ylab = "Expression Level (log2)",
        col = c("lightblue", "pink"),
        border = "black", notch = TRUE)

# 添加統計顯著性檢驗 (t檢驗)
urticaria <- data$Expression[data$Group == "Urticaria"]
normal_skin <- data$Expression[data$Group == "Normal Skin"]
p_value <- t.test(urticaria, normal_skin)$p.value
legend("topright", legend = paste("p-value =", signif(p_value, 3)), bty = "n")
