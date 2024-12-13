##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
if(!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if(!require('Seurat')) install.packages('Seurat'); library(Seurat)

#### Load data ####
Set_Dataset = "GSE202352_Normal"
if(Set_Dataset == "KGD_HS23378"){
  Seurat_Object <-  readRDS(file ="C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/KGL_Lab_EPPK_HS23378_2862.rds")
  Seurat_Object$Dataset <- "KGD_EPPK_HS23378_2862"
  
}else if(Set_Dataset == "GSE202352_Normal"){
  Seurat_Object <-  readRDS(file ="C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/GSE202352_combined_seurat_integrated.rds")
  Seurat_Object$Dataset <- "GSE202352_Normal"
  
}else if(Set_Dataset == "GSE189889_Acral_melanoma"){
  Seurat_Object <-  readRDS(file ="C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/GSE189889_Acral_Melanoma_SeuratObject_KGDLab.rds")
  Seurat_Object$Dataset <- "GSE189889_Acral_melanoma"
  
}else if(Set_Dataset == "GSE215121_Acral_melanoma"){
  Seurat_Object <-  readRDS(file ="C:/Users/q2330/Dropbox/##_GitHub/##_KGD_Lab/Acral_Melanoma/Export/GSE215121_Acral_Melanoma_SeuratObject_Integrated.rds")
  Seurat_Object$Dataset <- "GSE215121_Acral_melanoma"
  
}

DimPlot(Seurat_Object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

DimPlot(Seurat_Object, reduction = "umap",group.by = "Dataset", label = TRUE, pt.size = 0.5) + NoLegend()

Marker.lt <- c("EPGN", "EGFR", "KRT1", "KRT9", 
               "KIT", "EP300", "NF1", "COL1A1")

VlnPlot(Seurat_Object, features = Marker.lt, ncol = 4)

################################################################################
library(RColorBrewer)
# Generate a consistent color palette for unique clusters
color_palette <- colorRampPalette(brewer.pal(8, "Set3"))(length(unique(Seurat_Object$seurat_clusters)))

# DimPlot with consistent colors, black border, and increased text size
DimPlot(Seurat_Object, reduction = "umap", label = TRUE,
        group.by = "seurat_clusters", cols = color_palette) + 
  theme_minimal() +
  ggtitle("UMAP Plot with Seurat Clusters") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # Title size
    axis.title = element_text(size = 16, face = "bold"), # Axis title size
    axis.text = element_text(size = 14), # Axis text size
    legend.title = element_text(size = 16), # Legend title size
    legend.text = element_text(size = 14), # Legend text size
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5) # Black border around the plot
  )

# VlnPlot with consistent colors
DefaultAssay(Seurat_Object) <- "RNA"
VlnPlot(Seurat_Object, features = Marker.lt, ncol = 4, pt.size = 0.2, cols = color_palette) + 
  theme_minimal() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  # ggtitle("Expression of Selected Markers with Consistent Cluster Colors") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )



library(Seurat)
library(patchwork)

# Generate individual plots and store them in a list
plots <- lapply(Marker.lt, function(marker) {
  p <- VlnPlot(Seurat_Object, features = marker, pt.size = 0, cols = color_palette) + # No default points
    geom_jitter(aes(y = .data[[marker]]), width = 0.2, alpha = 0.1, size = 0.01) + # Transparent points for each plot
    theme_minimal() +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    ggtitle(paste( marker)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1) # Tilt x-axis labels by 45 degrees
    ) +
    NoLegend() # Remove the legend
  return(p)
})

# Combine plots into a single grid
wrap_plots(plots, ncol = 4)
