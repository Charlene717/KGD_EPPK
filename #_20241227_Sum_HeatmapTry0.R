#### Heatmap1 ####

# 確認並安裝必要套件
if(!require('ComplexHeatmap')) {install.packages('ComplexHeatmap'); library(ComplexHeatmap)}
if(!require('circlize')) {install.packages('circlize'); library(circlize)}
if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}
if(!require('tidyr')) {install.packages('tidyr'); library(tidyr)}
if(!require('tibble')) {install.packages('tibble'); library(tibble)}

# 處理資料：計算-log10(p_value)並標記顯著性
combined_data <- combined_data %>%
  mutate(
    log10_p_value = -log10(p_value),
    x_group = paste(Marker, Group1, sep = "-"),
    y_group = paste(Disease, Dataset, DataType, sep = "-"),
    is_significant = ifelse(p_value < 0.05, "*", "") # p_value < 0.05 則標記星號
  )

# 創建矩陣（Heatmap 使用矩陣格式）
heatmap_matrix <- combined_data %>%
  dplyr::select(x_group, y_group, log2FC) %>%
  tidyr::pivot_wider(names_from = x_group, values_from = log2FC) %>%
  tibble::column_to_rownames("y_group") %>%
  as.matrix()

# 創建星號矩陣
star_matrix <- combined_data %>%
  dplyr::select(x_group, y_group, is_significant) %>%
  tidyr::pivot_wider(names_from = x_group, values_from = is_significant) %>%
  tibble::column_to_rownames("y_group") %>%
  as.matrix()

# 替換 NA 值（填充為空字串）
star_matrix[is.na(star_matrix)] <- ""

# Heatmap 配色（包含灰色表示 NA）
heatmap_colors <- colorRamp2(
  c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  c("blue", "white", "red")
)

# 自訂背景色，用於填充 NA
na_color <- "grey"

# 繪製 Heatmap，關閉行列聚類
heatmap <- Heatmap(
  heatmap_matrix,
  name = "log2FC",
  col = heatmap_colors,
  na_col = na_color, # 指定 NA 的顏色為灰色
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  cluster_rows = FALSE, # 關閉行聚類
  cluster_columns = FALSE, # 關閉列聚類
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(star_matrix[i, j], x, y, gp = gpar(fontsize = 12, col = "black"))
  }
)

# 顯示 Heatmap
draw(heatmap)



#################################################################################
# 
# #### Heatmap2 ####
# 
# # 確認並安裝必要套件
# if(!require('ComplexHeatmap')) {install.packages('ComplexHeatmap'); library(ComplexHeatmap)}
# if(!require('circlize')) {install.packages('circlize'); library(circlize)}
# if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}
# if(!require('tidyr')) {install.packages('tidyr'); library(tidyr)}
# if(!require('tibble')) {install.packages('tibble'); library(tibble)}
# 
# # 處理資料：計算-log10(p_value)並標記顯著性
# combined_data <- combined_data %>%
#   mutate(
#     log10_p_value = -log10(p_value),
#     x_group = paste(Marker, Group1, sep = "-"),
#     y_group = paste(Disease, Dataset, DataType, sep = "-"),
#     is_significant = ifelse(p_value < 0.05, "*", "") # p_value < 0.05 則標記星號
#   )
# 
# # 創建矩陣（Heatmap 使用矩陣格式）
# heatmap_matrix <- combined_data %>%
#   dplyr::select(x_group, y_group, log2FC) %>%
#   tidyr::pivot_wider(names_from = x_group, values_from = log2FC) %>%
#   tibble::column_to_rownames("y_group") %>%
#   as.matrix()
# 
# # 替換 NA 值（填充為 0，或根據需要選擇其他值）
# heatmap_matrix[is.na(heatmap_matrix)] <- 0
# 
# # 創建星號矩陣
# star_matrix <- combined_data %>%
#   dplyr::select(x_group, y_group, is_significant) %>%
#   tidyr::pivot_wider(names_from = x_group, values_from = is_significant) %>%
#   tibble::column_to_rownames("y_group") %>%
#   as.matrix()
# 
# # 替換 NA 值（填充為空字串）
# star_matrix[is.na(star_matrix)] <- ""
# 
# # Heatmap 配色
# heatmap_colors <- colorRamp2(
#   c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
#   c("blue", "white", "red")
# )
# 
# # 繪製 Heatmap
# heatmap <- Heatmap(
#   heatmap_matrix,
#   name = "log2FC",
#   col = heatmap_colors,
#   row_names_gp = gpar(fontsize = 10, fontface = "bold"),
#   column_names_gp = gpar(fontsize = 10, fontface = "bold"),
#   show_row_dend = FALSE,
#   show_column_dend = FALSE,
#   cell_fun = function(j, i, x, y, width, height, fill) {
#     grid.text(star_matrix[i, j], x, y, gp = gpar(fontsize = 12, col = "black"))
#   }
# )
# 
# # 顯示 Heatmap
# draw(heatmap)
# 
