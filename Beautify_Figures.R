##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
if(!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if(!require('Seurat')) install.packages('Seurat'); library(Seurat)

#### Load data ####

KGD_EPPK_HS23378_2862_Seurat <-  readRDS(file ="C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/KGL_Lab_EPPK_HS23378_2862.rds")
DimPlot(KGD_EPPK_HS23378_2862_Seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
KGD_EPPK_HS23378_2862_Seurat$Dataset <- "KGD_EPPK_HS23378_2862"
DimPlot(KGD_EPPK_HS23378_2862_Seurat, reduction = "umap",group.by = "Dataset", label = TRUE, pt.size = 0.5) + NoLegend()

Marker.lt <- c("EPGN", "EGFR", "KRT1", "KRT9", 
               "KIT", "EP300", "NF1", "COL1A1")

VlnPlot(KGD_EPPK_HS23378_2862_Seurat, features = Marker.lt, ncol = 4)

################################################################################
# Generate a consistent color palette for unique clusters
color_palette <- colorRampPalette(brewer.pal(8, "Set3"))(length(unique(KGD_EPPK_HS23378_2862_Seurat$seurat_clusters)))

# DimPlot with consistent colors, black border, and increased text size
DimPlot(KGD_EPPK_HS23378_2862_Seurat, reduction = "umap", label = TRUE,
        group.by = "seurat_clusters", cols = color_palette) + 
  theme_minimal() +
  ggtitle("UMAP Plot with Consistent Cluster Colors") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # Title size
    axis.title = element_text(size = 16, face = "bold"), # Axis title size
    axis.text = element_text(size = 14), # Axis text size
    legend.title = element_text(size = 16), # Legend title size
    legend.text = element_text(size = 14), # Legend text size
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5) # Black border around the plot
  )

# VlnPlot with consistent colors
VlnPlot(KGD_EPPK_HS23378_2862_Seurat, features = Marker.lt, ncol = 4, pt.size = 0.2, cols = color_palette) + 
  theme_minimal() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  # ggtitle("Expression of Selected Markers with Consistent Cluster Colors") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

