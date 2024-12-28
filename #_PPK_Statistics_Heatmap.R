###############################################################################
#### 1. 載入 / 安裝套件 ####
###############################################################################
# 若未安裝，請透過 install.packages("xxx") 安裝
required_packages <- c('dplyr', 'readr', 'ggplot2', 'scales', 
                       'ComplexHeatmap', 'circlize', 'tidyr', 'tibble')

installed_packages <- rownames(installed.packages())

for(pkg in required_packages){
  if(!pkg %in% installed_packages){
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

###############################################################################
#### 2. 讀取檔案並整合資料 ####
###############################################################################
# 設定資料夾路徑 (請依實際情況修改)
data_folder <- "C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/KGD_EPPK/Export_Sum_20241228"

# 取得資料夾中所有 _stats.csv 檔案的路徑
file_list <- list.files(path = data_folder, pattern = "_stats\\.csv$", full.names = TRUE)

# 讀取並整合所有檔案
combined_data <- file_list %>%
  lapply(read_csv) %>%
  bind_rows()

# 重命名 Group1 為 Control
colnames(combined_data)[colnames(combined_data) == "Group1"] <- "Control"

# 檢視整合後的資料框（僅顯示前 6 列）
print(head(combined_data))

# 若需將整合後的資料另存為新檔案，可使用：
# output_file <- file.path(data_folder, "combined_stats.csv")
# write_csv(combined_data, output_file)

###############################################################################
#### 3. 繪製 Bubble Plot ####
###############################################################################
combined_data <- combined_data %>%
  mutate(
    log10_p_adjusted = -log10(p_adjusted), 
    x_group = Marker,  # 將 x_group 改為僅 Marker
    y_group = paste(Disease, Dataset, DataType, Control, sep = "-")  # 將 Control 加入 y_group
  )

p <- ggplot(combined_data, aes(x = x_group, y = y_group)) +
  geom_point(aes(size = log10_p_adjusted, color = log2FC), alpha = 0.7) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "log2FC"
  ) +
  scale_size_continuous(
    name = "-log10(p_adjusted)",
    range = c(2, 10)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid.major = element_line(color = "gray80", linetype = "dotted"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right"
  ) +
  labs(
    x = "Marker",
    y = "Disease, Dataset, DataType, and Control",
    title = "Bubble Plot of Marker Comparisons",
    subtitle = "Bubble size: -log10(p_adjusted); Color: log2FC"
  )

print(p)

###############################################################################
#### 4. 繪製 Heatmap（確保 Annotation 與主圖對齊）####
###############################################################################
# 1) 預處理：計算 -log10(p_value)，並標記顯著性
combined_data <- combined_data %>%
  mutate(
    log10_p_value = -log10(p_value),
    is_significant = ifelse(p_value < 0.05, "*", "")
  )

# 2) 建立 heatmap_matrix（行: y_group，列: x_group）
heatmap_matrix <- combined_data %>%
  dplyr::select(x_group, y_group, log2FC) %>%
  tidyr::pivot_wider(names_from = x_group, values_from = log2FC) %>%
  tibble::column_to_rownames("y_group") %>%
  as.matrix()

# 3) 建立星號矩陣 star_matrix
star_matrix <- combined_data %>%
  dplyr::select(x_group, y_group, is_significant) %>%
  tidyr::pivot_wider(names_from = x_group, values_from = is_significant) %>%
  tibble::column_to_rownames("y_group") %>%
  as.matrix()

# NA 補空字串
star_matrix[is.na(star_matrix)] <- ""

# -----------------------------------------------------------------------------
# 自定義顏色 (可依實際需求增減)
# -----------------------------------------------------------------------------
custom_colors <- list(
  Disease = c(
    "Atopic Dermatitis" = "#408080",
    "Lichen Planus"     = "#0f7344",
    "Urticaria"         = "#808040",
    "Psoriasis"         = "#9F4D95",
    "Ichthyosis_LI"     = "#f5d3d3", 
    "Ichthyosis_EI"     = "#fcb3b3", 
    "Ichthyosis_CIE"    = "#f78383", 
    "Ichthyosis_NS"     = "#e85a5a",
    "EPPK"              = "#BF0060", 
    "NPPK"              = "#930000"
  ),
  Dataset = c(
    "KGD_Lab"   = "#AE00AE",
    "GSE192832" = "#BE77FF",
    "GSE54456"  = "#D3A4FF",
    "GSE57178"  = "#3f88d1",
    "GSE206391" = "#005AB5"
  ),
  DataType = c(
    "RNA-seq"               = "#6C6C6C",
    "Array"                 = "#BEBEBE",
    "Spatial transcriptomics"= "#3C3C3C"
  ),
  Marker = c(
    "EPGN" = "#D9006C",
    "EGFR" = "#FF79BC"
  ),
  Control = c(
    "Non-lesional" = "#81C0C0",
    "Healthy"      = "#408080"
  )
)

# 可以自訂想要的 Disease 與 Dataset 排序
disease_order <- c(
  "Atopic Dermatitis", "Lichen Planus", "Urticaria", "Psoriasis",
  "Ichthyosis_LI", "Ichthyosis_EI", "Ichthyosis_CIE", "Ichthyosis_NS",
  "EPPK", "NPPK"
)
dataset_order <- c("KGD_Lab", "GSE192832", "GSE54456", "GSE57178", "GSE206391")

# -----------------------------------------------------------------------------
# 4) 建立 row_annotation_data，並且根據「自訂順序」排序
# -----------------------------------------------------------------------------
row_annotation_data <- combined_data %>%
  dplyr::select(y_group, Disease, Dataset, DataType, Control) %>%
  distinct() %>%
  # 轉為 factor 並設定順序
  mutate(
    Disease = factor(Disease, levels = disease_order),
    Dataset = factor(Dataset, levels = dataset_order),
    Control = factor(Control, levels = c("Non-lesional", "Healthy"))  # 設定 Control 的順序
  ) %>%
  arrange(Disease, Dataset, Control) %>%
  tibble::column_to_rownames("y_group")

# -----------------------------------------------------------------------------
# 5) 建立 column_annotation_data，並定義列順序 (Marker)
# -----------------------------------------------------------------------------
column_annotation_data <- combined_data %>%
  dplyr::select(x_group, Marker) %>%
  distinct() %>%
  tibble::column_to_rownames("x_group")

# -----------------------------------------------------------------------------
# 6) 以 row_annotation_data、column_annotation_data 中的 rownames() 為基準，
#    來「重新排序」heatmap_matrix、star_matrix，使 Heatmap 與 Annotation 一致
# -----------------------------------------------------------------------------
row_order <- rownames(row_annotation_data)
col_order <- rownames(column_annotation_data)

# 為避免 NA 或多餘行列，先檢查取交集 (通常是整齊對應時即可直接用 row_order, col_order)
row_order <- intersect(row_order, rownames(heatmap_matrix))
col_order <- intersect(col_order, colnames(heatmap_matrix))

# 重新對應排序
heatmap_matrix <- heatmap_matrix[row_order, col_order, drop = FALSE]
star_matrix   <- star_matrix[row_order, col_order, drop = FALSE]

# -----------------------------------------------------------------------------
# 7) 建立 rowAnnotation()、columnAnnotation()
# -----------------------------------------------------------------------------
row_annotation <- rowAnnotation(
  Disease = row_annotation_data$Disease,
  Dataset = row_annotation_data$Dataset,
  DataType = row_annotation_data$DataType,
  Control = row_annotation_data$Control,
  col = list(
    Disease = custom_colors$Disease,
    Dataset = custom_colors$Dataset,
    DataType = custom_colors$DataType,
    Control = custom_colors$Control  # 添加 Control 的顏色
  )
)

column_annotation <- columnAnnotation(
  Marker = column_annotation_data$Marker,
  col = list(
    Marker = custom_colors$Marker
  )
)

# -----------------------------------------------------------------------------
# 8) 設定 Heatmap 顏色階梯
# -----------------------------------------------------------------------------
# heatmap_colors <- colorRamp2(
#   c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
#   c("#465dcf", "white", "#e64343")
# )

if(!require('circlize')) {install.packages('circlize'); library(circlize)}

# 定義更多的分段點和對應的顏色
heatmap_colors <- colorRamp2(
  c(
    -2,
    # min(heatmap_matrix, na.rm = TRUE), 
    quantile(heatmap_matrix, 0.15, na.rm = TRUE),
    0, 
    quantile(heatmap_matrix, 0.85, na.rm = TRUE), 
    max(heatmap_matrix, na.rm = TRUE)
  ),
  c("#9dbee3", "#e1edfa", "white", "#f2c4c4", "#c22d2d")
)

# 使用 heatmap_colors 來生成熱圖

# -----------------------------------------------------------------------------
# 9) 繪製 Heatmap
# -----------------------------------------------------------------------------
heatmap <- Heatmap(
  heatmap_matrix,
  name = "log2FC",
  col = heatmap_colors,
  na_col = "grey", 
  show_row_names    = FALSE, 
  show_column_names = FALSE,
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    # 在每格繪製星號
    grid.text(star_matrix[i, j], x, y, gp = gpar(fontsize = 12, col = "black"))
  },
  top_annotation  = column_annotation,
  left_annotation = row_annotation
)

# 顯示 Heatmap
draw(heatmap)


heatmap2 <- Heatmap(
  heatmap_matrix,
  name = "log2FC",
  col = heatmap_colors,
  na_col = "grey", 
  show_row_names    = FALSE, 
  show_column_names = TRUE,
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    # 在每格繪製星號
    grid.text(star_matrix[i, j], x, y, gp = gpar(fontsize = 12, col = "black"))
  },
  column_names_side = "top",  # 將 x 軸文字移到上方
  column_names_gp = gpar(fontsize = 14, fontface = "bold"),  # 設定字型大小和粗體
  column_names_rot = 0,  # 設定文字方向為水平（橫向）
  column_title_rot = 0,  # 取消標題的旋轉
  column_names_centered = TRUE,  # 文字相對於自己的欄位置中
  left_annotation = row_annotation  # 移除 top_annotation
)


heatmap2

###############################################################################
#### 5. 輸出檔案 ####
###############################################################################
# 建立輸出資料夾與檔名
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 10) 
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))
Set_note <- paste0(Name_FileID, "_Sum")
Name_Export <- paste0(Name_FileID)
Name_ExportFolder <- paste0("Export_", Set_note)

if (!dir.exists(Name_ExportFolder)){
  dir.create(Name_ExportFolder)
}

# 輸出 PDF
pdf(file = paste0(Name_ExportFolder, "/", Name_Export, "_heatmap.pdf"),
    width = 8, height = 8)
print(heatmap2)
print(heatmap)
dev.off()

# 輸出整個 R 環境 (包含 combined_data、圖形物件等)
save.image(file = paste0(Name_ExportFolder, "/", Name_Export, ".RData"))
