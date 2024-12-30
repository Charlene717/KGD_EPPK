# 確保沒有負值，將負值設為零
combined_seurat <- ScaleData(combined_seurat, do.center = FALSE)  # 避免居中產生負值

# 繪製小提琴圖並限制最小值為零
plots <- VlnPlot(
  combined_seurat,
  features = c("EPGN", "EGFR"),
  split.by = "orig.ident",
  group.by = "seurat_clusters",
  pt.size = 0,
  combine = FALSE
)

# 自定義去除負值的方法：限制 y 軸最小值
plots <- lapply(plots, function(p) {
  p + ylim(0, NA)  # 將 y 軸的最小值設為 0，最大值自動計算
})

# 排列圖形
wrap_plots(plots = plots, ncol = 1)


######################################################################################
# 提取數據
library(ggplot2)
library(dplyr)
library(tidyr)

# 提取要繪製的數據
data_to_plot <- FetchData(combined_seurat, vars = c("EPGN", "EGFR", "seurat_clusters", "orig.ident")) %>%
  pivot_longer(
    cols = c("EPGN", "EGFR"), 
    names_to = "Feature", 
    values_to = "Expression"
  )

# 繪製小提琴圖
ggplot(data_to_plot, aes(x = seurat_clusters, y = Expression, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(~Feature, ncol = 1) +
  theme_minimal() +
  labs(x = "Clusters", y = "Expression Level", fill = "Origin") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
