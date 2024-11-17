##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
if(!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if(!require('Seurat')) install.packages('Seurat'); library(Seurat)
if(!require('hdf5r')) install.packages("hdf5r"); library(hdf5r)
if(!require('dplyr')) install.packages("dplyr"); library(dplyr)
if(!require('Matrix')) install.packages("Matrix"); library(Matrix)

#### Load Data ####
## Update the folder path to reflect new location
folder_path <- "D:/Dropbox/KGD_Lab/20241015 (待整理)(線上資料整理)_Acral Melanoma/Acral melanoma datasets/GSE215121_RAW/"

# Load all .h5 files in the directory
data_files <- list.files(path = folder_path, pattern = "\.h5$", full.names = TRUE)

# Initialize an empty list to store expression matrices and metadata
expression_matrices <- list()
metadata_list <- list()

# Loop through each .h5 file and extract the data
for (file in data_files) {
  # Extract file name information for metadata
  file_name <- basename(file)
  file_info <- strsplit(file_name, "_filtered")[[1]][1]
  file_metadata <- strsplit(file_info, "_")[[1]]
  sample_id <- file_metadata[1]
  condition <- file_metadata[2]
  
  # Store metadata information
  metadata <- data.frame(SampleID = sample_id, Condition = condition, stringsAsFactors = FALSE)
  metadata_list[[file]] <- metadata
  
  data <- H5File$new(file, mode = "r")
  data_matrix_group <- data[['matrix']]
  
  # ## 檢查用
  # # 列出 matrix group 中的所有項目
  # items_in_group <- data_matrix_group$ls()
  # 
  # # 獲取 features group
  # features_group <- data_matrix_group[['features']]
  # 
  # # 列出 features group 中的所有項目
  # items_in_features <- features_group$ls()
  
  # Read sparse matrix data
  indices <- data_matrix_group[['indices']]$read()
  indptr <- data_matrix_group[['indptr']]$read()
  data_values <- data_matrix_group[['data']]$read()
  matrix_shape <- data_matrix_group[['shape']]$read()
  
  # Reconstruct sparse matrix
  expr_matrix <- sparseMatrix(i = indices + 1, p = indptr, x = data_values, dims = c(matrix_shape[1], matrix_shape[2]))
  
  # Add row and column names if they exist and match the dimensions
  features_group <- data_matrix_group[['features']]
  
  # Extract gene names
  gene_names <- tryCatch(features_group[['name']]$read(), error = function(e) NULL)
  cell_names <- tryCatch(data_matrix_group[['barcodes']]$read(), error = function(e) NULL)
  
  if (!is.null(gene_names) && length(gene_names) == nrow(expr_matrix)) {
    rownames(expr_matrix) <- gene_names
  } else {
    rownames(expr_matrix) <- paste0("Gene", seq_len(nrow(expr_matrix))) # Fallback gene names
    warning("Gene names are missing or do not match matrix dimensions. Using placeholder names.")
  }
  
  if (!is.null(cell_names) && length(cell_names) == ncol(expr_matrix)) {
    colnames(expr_matrix) <- cell_names
  } else {
    colnames(expr_matrix) <- paste0("Cell", seq_len(ncol(expr_matrix))) # Fallback cell names
    warning("Cell names are missing or do not match matrix dimensions. Using placeholder names.")
  }
  
  # Append the expression matrix to the list
  expression_matrices[[file]] <- expr_matrix
  
  data$close_all() # Close the file to free resources
}

# Merge all matrices into a single matrix
combined_matrix <- do.call(cbind, expression_matrices)

# Merge all metadata into a single data frame if available
if (length(metadata_list) > 0) {
  combined_metadata <- do.call(rbind, metadata_list)
} else {
  combined_metadata <- NULL
  message("No metadata available to merge.")
}

# Check for duplicated cell names
duplicate_cells <- duplicated(colnames(combined_matrix))
if (any(duplicate_cells)) {
  message("Duplicated cell names found.")
  print(colnames(combined_matrix)[duplicate_cells])
} else {
  message("No duplicated cell names found.")
}

# Ensure unique gene and cell names
rownames(combined_matrix) <- make.unique(rownames(combined_matrix))
colnames(combined_matrix) <- make.unique(colnames(combined_matrix))

# Ensure no NA values in gene and cell names
if (any(is.na(rownames(combined_matrix)))) {
  message("NA values found in gene names, replacing with unique names.")
  rownames(combined_matrix)[is.na(rownames(combined_matrix))] <- paste0("Gene", seq_len(sum(is.na(rownames(combined_matrix)))))
}

if (any(is.na(colnames(combined_matrix)))) {
  message("NA values found in cell names, replacing with unique names.")
  colnames(combined_matrix)[is.na(colnames(combined_matrix))] <- paste0("Cell", seq_len(sum(is.na(colnames(combined_matrix)))))
}

#### Create Seurat Object ####
acral_melanoma <- CreateSeuratObject(counts = combined_matrix, project = "AcralMelanoma", min.cells = 3, min.features = 200)

# Add metadata to Seurat object if available
if (!is.null(combined_metadata)) {
  acral_melanoma <- AddMetaData(acral_melanoma, metadata = combined_metadata)
}

# Add mitochondrial percentage to metadata
acral_melanoma[["percent.mt"]] <- PercentageFeatureSet(acral_melanoma, pattern = "^MT-")

# Quality control filtering
acral_melanoma <- subset(acral_melanoma, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Data normalization
acral_melanoma <- NormalizeData(acral_melanoma)

# Identify highly variable features
acral_melanoma <- FindVariableFeatures(acral_melanoma, selection.method = "vst", nfeatures = 2000)

# Scale the data
all.genes <- rownames(acral_melanoma)
acral_melanoma <- ScaleData(acral_melanoma, features = all.genes)

# Principal Component Analysis (PCA)
acral_melanoma <- RunPCA(acral_melanoma, features = VariableFeatures(object = acral_melanoma))

# Visualize PCA results
ElbowPlot(acral_melanoma)

# Clustering analysis using PCA
acral_melanoma <- FindNeighbors(acral_melanoma, dims = 1:10)
acral_melanoma <- FindClusters(acral_melanoma, resolution = 0.5)

# UMAP for dimensionality reduction
acral_melanoma <- RunUMAP(acral_melanoma, dims = 1:10)

# Visualize UMAP results
DimPlot(acral_melanoma, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(acral_melanoma, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "curated.cell.types") + NoLegend()
DimPlot(acral_melanoma, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "sample") + NoLegend()

# Find markers for each cluster
cluster_markers <- FindAllMarkers(acral_melanoma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Display marker genes
print(cluster_markers)

DefaultAssay(acral_melanoma) <- "RNA"

Marker.lt <- c("EPGN", "EGFR","KRT1","KRT16",
               "BRAF","NRAS","NF1","TERT")

Marker.lt <- c("EPGN", "EGFR","KRT1","KRT9",
               "KIT","EP300","NF1","COL1A1")
VlnPlot(acral_melanoma, features = Marker.lt,ncol = 4)
VlnPlot(acral_melanoma, features = Marker.lt[1:4],ncol = 2, group.by = "blueprint.main")
VlnPlot(acral_melanoma, features = Marker.lt[5:8],ncol = 2, group.by = "blueprint.main")
FeaturePlot(acral_melanoma, features = Marker.lt,ncol = 4)

FeaturePlot(acral_melanoma, features = c("EPGN","EGFR","NPEPL1","SERPINA12"),ncol = 4)

# Save Seurat object
saveRDS(acral_melanoma, file = paste0(folder_path, "GSE215121_Acral_Melanoma_SeuratObject.rds"))
