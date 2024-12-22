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
folder_path <- "C:/Charlene/Dataset_Online/GSE206391/"

# Load all .h5 files in the directory
data_files <- list.files(path = folder_path, pattern = "\\.h5$", full.names = TRUE)

# Initialize lists to store Seurat objects and metadata
seurat_objects <- list()
metadata_list <- list()

# Loop through each .h5 file, create Seurat object, and extract metadata
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
    rownames(expr_matrix) <- make.unique(gene_names)
  } else {
    rownames(expr_matrix) <- paste0("Gene", seq_len(nrow(expr_matrix))) # Fallback gene names
    warning("Gene names are missing or do not match matrix dimensions. Using placeholder names.")
  }
  
  if (!is.null(cell_names) && length(cell_names) == ncol(expr_matrix)) {
    colnames(expr_matrix) <- make.unique(cell_names)
  } else {
    colnames(expr_matrix) <- paste0("Cell", seq_len(ncol(expr_matrix))) # Fallback cell names
    warning("Cell names are missing or do not match matrix dimensions. Using placeholder names.")
  }
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr_matrix, project = sample_id, min.cells = 3, min.features = 200)
  
  # Add metadata to Seurat object
  metadata <- metadata[rep(1, ncol(expr_matrix)), ]
  rownames(metadata) <- colnames(expr_matrix)
  seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)
  
  # Add mitochondrial percentage to metadata
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Quality control filtering
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
  
  # Append Seurat object to the list
  seurat_objects[[file]] <- seurat_obj
  
  data$close_all() # Close the file to free resources
}

# Normalize and find variable features for each Seurat object
for (i in 1:length(seurat_objects)) {
  seurat_objects[[i]] <- NormalizeData(seurat_objects[[i]])
  seurat_objects[[i]] <- FindVariableFeatures(seurat_objects[[i]], selection.method = "vst", nfeatures = 2000)
}

# Integrate Seurat objects
anchors <- FindIntegrationAnchors(object.list = seurat_objects, dims = 1:20)
acral_melanoma_integrated <- IntegrateData(anchorset = anchors, dims = 1:20)

# Scale the data
all.genes <- rownames(acral_melanoma_integrated)
acral_melanoma_integrated <- ScaleData(acral_melanoma_integrated, features = all.genes)

# Principal Component Analysis (PCA)
acral_melanoma_integrated <- RunPCA(acral_melanoma_integrated, features = VariableFeatures(object = acral_melanoma_integrated))

# Visualize PCA results
ElbowPlot(acral_melanoma_integrated)

# Clustering analysis using PCA
acral_melanoma_integrated <- FindNeighbors(acral_melanoma_integrated, dims = 1:10)
acral_melanoma_integrated <- FindClusters(acral_melanoma_integrated, resolution = 0.5)

# UMAP for dimensionality reduction
acral_melanoma_integrated <- RunUMAP(acral_melanoma_integrated, dims = 1:10)

# Visualize UMAP results
DimPlot(acral_melanoma_integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(acral_melanoma_integrated, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "curated.cell.types") + NoLegend()
DimPlot(acral_melanoma_integrated, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "sample") + NoLegend()

# Find markers for each cluster
cluster_markers <- FindAllMarkers(acral_melanoma_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Display marker genes
print(cluster_markers)

DefaultAssay(acral_melanoma_integrated) <- "RNA"

Marker.lt <- c("EPGN", "EGFR","KRT1","KRT16",
               "BRAF","NRAS","NF1","TERT")

Marker.lt <- c("EPGN", "EGFR","KRT1","KRT9",
               "KIT","EP300","NF1","COL1A1")
VlnPlot(acral_melanoma_integrated, features = Marker.lt,ncol = 4)
VlnPlot(acral_melanoma_integrated, features = Marker.lt[1:4],ncol = 2, group.by = "blueprint.main")
VlnPlot(acral_melanoma_integrated, features = Marker.lt[5:8],ncol = 2, group.by = "blueprint.main")
FeaturePlot(acral_melanoma_integrated, features = Marker.lt,ncol = 4)

FeaturePlot(acral_melanoma_integrated, features = c("EPGN","EGFR","NPEPL1","SERPINA12"),ncol = 4)

# Save Seurat object
saveRDS(acral_melanoma_integrated, file = paste0(folder_path, "GSE215121_Acral_Melanoma_SeuratObject_Integrated.rds"))
