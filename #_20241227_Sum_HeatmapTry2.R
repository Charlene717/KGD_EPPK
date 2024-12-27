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

# 動態生成顏色映射，涵蓋所有 Disease 的值
unique_diseases <- unique(combined_data$Disease)
custom_colors <- list(
  # Disease = setNames(scales::hue_pal()(length(unique_diseases)), unique_diseases),
  Disease = c("Atopic Dermatitis" = "#408080", "Lichen Planus" = "#0066CC", "Urticaria" = "#808040", 
              "Psoriasis" = "#9F4D95",
              "Ichthyosis_LI" = "#f5d3d3", "Ichthyosis_EI" = "#fcb3b3", "Ichthyosis_CIE" = "#f78383", "Ichthyosis_NS" = "#e85a5a",
              "EPPK" = "#BF0060", "NPPK" = "#D94600"),
  
  Dataset = c("GSE192832" = "#BE77FF", "GSE206391" = "#0066CC", "GSE54456" = "#FF0080",
              "GSE57178" = "#019858", "KGD_Lab" = "#BF0060"),
  DataType = c("RNA-seq" = "#5B5B5B", "Array" = "#9D9D9D", "Spatial transcriptomics" = "#3C3C3C"),
  Marker = c("EPGN" = "#FF79BC", "EGFR" = "#D9006C"),
  Group1 = c("Non-lesional" = "#02DF82", "Healthy" = "#009393")
)

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
    Disease = custom_colors$Disease,
    Dataset = custom_colors$Dataset,
    DataType = custom_colors$DataType
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
    Marker = custom_colors$Marker,
    Group1 = custom_colors$Group1
  )
)

# Heatmap 配色（包含灰色表示 NA）
heatmap_colors <- colorRamp2(
  c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  c("blue", "white", "red")
)

# 繪製 Heatmap
heatmap <- Heatmap(
  heatmap_matrix,
  name = "log2FC",
  col = heatmap_colors,
  na_col = "grey", # 指定 NA 的顏色為灰色
  show_row_names = FALSE, # 關閉行文字
  show_column_names = FALSE, # 關閉列文字
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
