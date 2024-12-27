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

# 行標註數據
row_annotation_data <- combined_data %>%
  dplyr::select(y_group, Disease, Dataset, DataType) %>%
  distinct() %>%
  tibble::column_to_rownames("y_group")

# 創建行的 Annotation
row_annotation <- rowAnnotation(
  Disease = row_annotation_data$Disease,
  Dataset = row_annotation_data$Dataset,
  DataType = row_annotation_data$DataType,
  col = list(
    Disease = structure(scales::hue_pal()(length(unique(row_annotation_data$Disease))),
                        names = unique(row_annotation_data$Disease)),
    Dataset = structure(scales::hue_pal()(length(unique(row_annotation_data$Dataset))),
                        names = unique(row_annotation_data$Dataset)),
    DataType = structure(scales::hue_pal()(length(unique(row_annotation_data$DataType))),
                         names = unique(row_annotation_data$DataType))
  )
)

# 列標註數據
column_annotation_data <- combined_data %>%
  dplyr::select(x_group, Marker, Group1) %>%
  distinct() %>%
  tibble::column_to_rownames("x_group")

# 創建列的 Annotation
column_annotation <- columnAnnotation(
  Marker = column_annotation_data$Marker,
  Group1 = column_annotation_data$Group1,
  col = list(
    Marker = structure(scales::hue_pal()(length(unique(column_annotation_data$Marker))),
                       names = unique(column_annotation_data$Marker)),
    Group1 = structure(scales::hue_pal()(length(unique(column_annotation_data$Group1))),
                       names = unique(column_annotation_data$Group1))
  )
)

# 繪製 Heatmap
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
  },
  top_annotation = column_annotation, # 添加列的 Annotation
  left_annotation = row_annotation    # 添加行的 Annotation
)

# 顯示 Heatmap
draw(heatmap)
