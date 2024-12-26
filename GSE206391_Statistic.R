# 檢查並載入必要套件
if (!require('dplyr')) { install.packages('dplyr'); library(dplyr) }

# 提取數據
Marker.lt <- c("EPGN", "EGFR")
plot_data <- FetchData(seurat_object, vars = c(Marker.lt, "group", "groupLN"))

# 定義比較組別
group_comparisons <- list(
  c( "AD_NON LESIONAL", "AD_LESIONAL"),
  c("LP_NON LESIONAL", "LP_LESIONAL"),
  c("Pso_NON LESIONAL", "Pso_LESIONAL")
)

# 計算 logFC 和 p 值的函數
calculate_stats <- function(data, marker, groups) {
  # 篩選數據
  data_filtered <- data %>%
    filter(group %in% groups) %>%
    select(group, !!sym(marker))
  
  # 計算 logFC
  group_means <- data_filtered %>%
    group_by(group) %>%
    summarize(mean_expression = mean(!!sym(marker)), .groups = 'drop') %>%
    arrange(match(group, groups))
  
  logFC <- log2(group_means$mean_expression[2] / group_means$mean_expression[1])
  
  # 計算 p 值 (t 檢驗)
  group1_values <- data_filtered %>% filter(group == groups[1]) %>% pull(!!sym(marker))
  group2_values <- data_filtered %>% filter(group == groups[2]) %>% pull(!!sym(marker))
  p_value <- t.test(group1_values, group2_values)$p.value
  # p_value <- wilcox.test(group1_values, group2_values)$p.value
  
  
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
