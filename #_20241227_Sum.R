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
