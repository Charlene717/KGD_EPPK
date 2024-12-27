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
         x_group = paste(Marker, Group1, sep = "-"),
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
    x = "Marker and Group1",
    y = "Disease, Dataset, and DataType",
    title = "Bubble Plot of Marker Comparisons",
    subtitle = "Bubble size: -log10(p_adjusted); Color: log2FC"
  )

# 顯示圖表
print(p)

