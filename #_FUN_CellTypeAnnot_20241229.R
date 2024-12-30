##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("devtools")) install.packages("devtools"); library(devtools)
# if (!requireNamespace("SeuratData", quietly = TRUE)) {install.packages("SeuratData")}; library(SeuratData)
if (!requireNamespace("SeuratData", quietly = TRUE)) {devtools::install_github('satijalab/seurat-data')}; library(SeuratData)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("SingleCellExperiment"))  BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
if(!require("scater")) BiocManager::install("scater") ; library(scater)
if(!require("celldex")) BiocManager::install("celldex"); library(celldex)
if(!require("scuttle")) BiocManager::install("scuttle"); library(scuttle)


#### singleR ####
Run_singleR <- function(Query_Seurat, Reference_Seurat,
                        Set_RefAnnoCol = "Actual_Cell_Type", seurat_version = "V5", ...) {
  # Load necessary packages
  # Load necessary packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!require("SingleR", quietly = TRUE)) BiocManager::install("SingleR"); library(SingleR)
  if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
  if (!require("scater", quietly = TRUE)) BiocManager::install("scater"); library(scater)
  
  
  # Handle Query_Seurat layers
  layers <- Layers(Query_Seurat, assay = "RNA")
  data_layers <- grep("^data\\.", layers, value = TRUE)
  
  if (length(data_layers) < 2) stop("Insufficient data layers found in Query_Seurat")
  
  # Extract and unify cell names
  data_list <- lapply(data_layers, function(layer) {
    data <- GetAssayData(Query_Seurat, assay = "RNA", layer = layer)
    colnames(data) <- sub("^.*_", "", colnames(data))
    return(data)
  })
  
  # Find union of genes and align data
  all_genes <- Reduce(union, lapply(data_list, rownames))
  aligned_data <- lapply(data_list, function(data) {
    full_matrix <- matrix(0, nrow = length(all_genes), ncol = ncol(data),
                          dimnames = list(all_genes, colnames(data)))
    full_matrix[rownames(data), ] <- as.matrix(data)
    return(full_matrix)
  })
  
  # Merge data and create SingleCellExperiment for query
  merged_data <- do.call(cbind, aligned_data)
  query_sce <- SingleCellExperiment(assays = list(counts = merged_data))
  
  # Preprocess reference dataset
  ref_sce <- as.SingleCellExperiment(Reference_Seurat)
  ref_sce$label <- ref_sce[[Set_RefAnnoCol]]
  ref_sce <- ref_sce[, !is.na(ref_sce$label)]
  ref_sce <- tryCatch({
    scater::logNormCounts(ref_sce)
  }, error = function(e) ref_sce)
  
  # Log-normalize query dataset
  query_sce <- tryCatch({
    scater::logNormCounts(query_sce[, colSums(counts(query_sce)) > 0])
  }, error = function(e) query_sce)
  
  # Run SingleR
  SingleR.lt <- SingleR(
    test = query_sce,
    ref = ref_sce,
    assay.type.test = 1,
    labels = ref_sce$label,
    de.method = "classic"
  )
  
  # Add SingleR results to Seurat metadata
  Query_Seurat@meta.data$label_singleR_NoReject <- SingleR.lt$labels
  Query_Seurat@meta.data[[paste0("label_singleR")]] <- ifelse(
    is.na(SingleR.lt$pruned.labels), "Unassign", SingleR.lt$pruned.labels
  )
  
  # Store additional SingleR metrics in Seurat object
  Query_Seurat@misc$CTAnnot <- list(
    singleR_Scores = SingleR.lt@listData[["scores"]],
    singleR_Delta = SingleR.lt@listData[["delta.next"]]
  )
  
  return(Query_Seurat)
}
