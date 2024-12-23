# 確保載入必要套件
if (!require('dplyr')) install.packages('dplyr')
library(dplyr)

# 定義基因列表和比較組別
Marker.lt <- c("EPGN","EGFR")  # 目前只有一個基因
group_comparisons <- list(
  c("Lesional", "Non-lesional"),
  c("Lesional", "Normal")
)

# 提取數據
plot_data <- data  # 使用已經整理好的 data 數據框
# 確認格式：plot_data 應該至少包含列 `Expression` 和 `Group`

# 計算 logFC 和 p 值的函數
calculate_stats <- function(data, marker, groups) {
  # 篩選數據
  data_filtered <- data %>%
    filter(Group %in% groups) %>%
    select(Group, Expression)
  
  # 計算 logFC
  group_means <- data_filtered %>%
    group_by(Group) %>%
    summarize(mean_expression = mean(Expression), .groups = 'drop') %>%
    arrange(match(Group, groups))
  
  logFC <- log2(group_means$mean_expression[2] / group_means$mean_expression[1])
  
  # 計算 p 值 (t 檢驗)
  group1_values <- data_filtered %>% filter(Group == groups[1]) %>% pull(Expression)
  group2_values <- data_filtered %>% filter(Group == groups[2]) %>% pull(Expression)
  p_value <- t.test(group1_values, group2_values)$p.value
  
  return(data.frame(Marker = marker, Group1 = groups[1], Group2 = groups[2], logFC = logFC, p_value = p_value))
}

# 計算每個 Marker 的 logFC 和 p 值
stats_results <- lapply(Marker.lt, function(marker) {
  do.call(rbind, lapply(group_comparisons, function(groups) {
    calculate_stats(plot_data, marker, groups)
  }))
}) %>%
  do.call(rbind, .)

# 添加校正後的 p 值
stats_results$p_adjusted <- p.adjust(stats_results$p_value, method = "BH")  # Benjamini-Hochberg 方法

# 查看結果
print(stats_results)

# 將結果保存到文件（可選）
write.csv(stats_results, "group_comparisons_stats_with_adjusted_p.csv", row.names = FALSE)
