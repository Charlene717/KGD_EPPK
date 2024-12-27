# 設定資料夾路徑
data_folder <- "C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/KGD_EPPK/Export_Sum_20241227"

# 檢查是否已安裝必要的套件，未安裝則進行安裝
if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}
if(!require('readr')) {install.packages('readr'); library(readr)}

# 取得資料夾中所有 _stats.csv 檔案的路徑
file_list <- list.files(path = data_folder, pattern = "_stats\\.csv$", full.names = TRUE)

# 讀取並整合所有檔案
combined_data <- file_list %>%
  lapply(read_csv) %>%    # 使用 readr::read_csv 讀取每個 CSV 檔案
  bind_rows()             # 將所有資料框合併為一個

# 重命名 Group1 為 Control
colnames(combined_data)[colnames(combined_data) == "Group1"] <- "Control"

# 檢視整合後的資料框
print(head(combined_data))

# # 若需將整合後的資料另存為新檔案，可使用以下程式碼
# output_file <- file.path(data_folder, "combined_stats.csv")
# write_csv(combined_data, output_file)


###############################################################################

# 檢查必要套件並安裝
if(!require('ggplot2')) {install.packages('ggplot2'); library(ggplot2)}
if(!require('scales')) {install.packages('scales'); library(scales)}
if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}

# 計算額外欄位
combined_data <- combined_data %>%
  mutate(log10_p_adjusted = -log10(p_adjusted), 
         x_group = paste(Marker, Control, sep = "-"),
         y_group = paste(Disease, Dataset, DataType, sep = "-"))

# 繪製氣泡圖
p <- ggplot(combined_data, aes(x = x_group, y = y_group)) +
  geom_point(aes(size = log10_p_adjusted, color = log2FC), alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "log2FC") +
  scale_size_continuous(name = "-log10(p_adjusted)", range = c(2, 10)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid.major = element_line(color = "gray80", linetype = "dotted"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right"
  ) +
  labs(
    x = "Marker and Control",
    y = "Disease, Dataset, and DataType",
    title = "Bubble Plot of Marker Comparisons",
    subtitle = "Bubble size: -log10(p_adjusted); Color: log2FC"
  )

# 顯示圖表
print(p)

################################################################################
#### Heatmap ####
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
    x_group = paste(Marker, Control, sep = "-"),
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

# 自定義顏色
custom_colors <- list(
  Disease = c("Atopic Dermatitis" = "#408080", "Lichen Planus" = "#0f7344", "Urticaria" = "#808040", 
              "Psoriasis" = "#9F4D95",
              "Ichthyosis_LI" = "#f5d3d3", "Ichthyosis_EI" = "#fcb3b3", "Ichthyosis_CIE" = "#f78383", "Ichthyosis_NS" = "#e85a5a",
              "EPPK" = "#BF0060", "NPPK" = "#930000"),
  Dataset = c("KGD_Lab" = "#AE00AE", "GSE192832" = "#BE77FF", "GSE54456" = "#D3A4FF",
              "GSE57178" = "#3f88d1","GSE206391" = "#005AB5" ),
  DataType = c("RNA-seq" = "#5B5B5B", "Array" = "#9D9D9D", "Spatial transcriptomics" = "#3C3C3C"),
  Marker = c("EPGN" = "#D9006C", "EGFR" = "#FF79BC"),
  Control = c("Non-lesional" = "#81C0C0", "Healthy" = "#408080")
)

# 自定義 Dataset 和 Disease 的順序
dataset_order <- c("KGD_Lab", "GSE192832", "GSE54456", "GSE57178", "GSE206391")
disease_order <- c("Atopic Dermatitis", "Lichen Planus", "Urticaria", "Psoriasis",
                   "Ichthyosis_LI", "Ichthyosis_EI", "Ichthyosis_CIE", "Ichthyosis_NS",
                   "EPPK", "NPPK")

# 行標註數據
row_annotation_data <- combined_data %>%
  dplyr::select(y_group, Disease, Dataset, DataType) %>%
  distinct() %>%
  mutate(
    Disease = factor(Disease, levels = disease_order), # 設定 Disease 順序
    Dataset = factor(Dataset, levels = dataset_order)  # 設定 Dataset 順序
  ) %>%
  arrange(Disease, Dataset) %>% # 根據自定義順序排序
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
  dplyr::select(x_group, Marker, Control) %>%
  distinct() %>%
  tibble::column_to_rownames("x_group")

# 創建列的 Annotation
column_annotation <- columnAnnotation(
  Marker = column_annotation_data$Marker,
  Control = column_annotation_data$Control,
  col = list(
    Marker = custom_colors$Marker,
    Control = custom_colors$Control
  )
)

# Heatmap 配色（包含灰色表示 NA）
heatmap_colors <- colorRamp2(
  c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  c("#465dcf", "white", "#e64343")
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


#### Export ####
# Set export parameters
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 10) # Generate unique time-based ID
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))
Set_note <- paste0(Name_FileID, "_Sum")
Name_Export <- paste0(Name_FileID)
Name_ExportFolder <- paste0("Export_", Set_note)
# Create export folder if it does not exist
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}

# 將圖形輸出到 PDF
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_heatmap.pdf"), width = 10, height = 8)
print(heatmap)
dev.off()

#### Export ####
## Export RData
save.image(paste0(Name_ExportFolder, "/", Name_Export, ".RData"))

