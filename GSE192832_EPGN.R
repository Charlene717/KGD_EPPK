# 確保必要套件已安裝
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)

# 讀取表達數據和註解數據
expression_data <- read.csv("C:/Users/q2330/Dropbox/KGD_Lab/20241129_(PPK) Bulk RNA分析/GSE192832/GSE192832_expression_v2.csv", row.names = 1)
annotation_data <- readxl::read_excel("C:/Users/q2330/Dropbox/KGD_Lab/20241129_(PPK) Bulk RNA分析/GSE192832/GSE192832_Annotation.xlsx")

# 篩選出 LS 和 Control 組樣本
ls_samples <- annotation_data %>%
  filter(`!Sample_characteristics_ch1...3` == "lesion: LS") %>%
  pull(ID)

control_samples <- annotation_data %>%
  filter(`!Sample_characteristics_ch1...3` == "lesion: Control") %>%
  pull(ID)

# 將樣本 ID 格式調整為表達數據的列名格式
ls_samples <- paste0("X", ls_samples)
control_samples <- paste0("X", control_samples)

# 提取 EPGN 表達數據
epgn_expression <- expression_data["EPGN", ]

# 整理 LS 和 Control 組的 EPGN 表達值
ls_expression <- epgn_expression[ls_samples]
control_expression <- epgn_expression[control_samples]

# 統計分析 (t-test)
t_test_result <- t.test(as.numeric(ls_expression), as.numeric(control_expression), alternative = "two.sided")

# 可視化表達數據
data_to_plot <- data.frame(
  Expression = c(as.numeric(ls_expression), as.numeric(control_expression)),
  Group = c(rep("LS", length(ls_expression)), rep("Control", length(control_expression)))
)

ggplot(data_to_plot, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  labs(
    title = "EPGN Expression in LS vs Control",
    subtitle = paste0("t-test p-value: ", signif(t_test_result$p.value, 3)),
    x = "Group",
    y = "EPGN Expression"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("LS" = "#FF9999", "Control" = "#9999FF"))




#### 視覺化與統計分析 ####

# 計算樣本數
ls_count <- length(ls_expression)
hs_count <- length(control_expression)

# 更新數據框，調整組別名稱和順序
data_to_plot <- data.frame(
  Expression = c(as.numeric(ls_expression), as.numeric(control_expression)),
  Group = factor(c(rep("LS", length(ls_expression)), rep("NS", length(control_expression))),
                 levels = c("LS", "NS")) # LS 在左，NS 在右
)

# 更新圖表
ggplot(data_to_plot, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_signif(
    comparisons = list(c("LS", "NS")), # 比較的組別
    map_signif_level = TRUE, # 顯示顯著性標註
    y_position = max(data_to_plot$Expression) + 5, # 框線位置
    color = "black", # 框線顏色
    textsize = 6, # 顯著性文字大小
    face = "bold",
    test = "t.test" # 手動指定使用 t 檢定
  ) +
  labs(
    title = "EPGN Expression in LS vs NS (GSE192832)", # 包含 GSE 編號
    x = "Group",
    y = "EPGN Expression",
    caption = paste0("LS: Ichthyosis (n = ", ls_count, "); NS: Normal Skin (n = ", hs_count, ").") # 詳細註解
  ) +
  scale_fill_manual(values = c("LS" = "#FF9999", "NS" = "#9999FF")) + # 設置顏色
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"), # X 軸文字加大加粗
    axis.title.y = element_text(size = 16, face = "bold"), # Y 軸文字加大加粗
    axis.text.x = element_text(size = 14, face = "bold"),  # X 軸標籤加大加粗
    axis.text.y = element_text(size = 14, face = "bold"),  # Y 軸標籤加大加粗
    plot.caption = element_text(hjust = 0, size = 10),
    legend.position = "none", # 移除圖例
    aspect.ratio = 1, # 設置長寬比為 1:1
    panel.border = element_rect(color = "black", fill = NA, size = 1) # 黑色框線對齊 X 和 Y 軸
  )

