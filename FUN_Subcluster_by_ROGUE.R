## Ref: https://htmlpreview.github.io/?https://github.com/PaulingLiu/ROGUE/blob/master/vignettes/ROGUE_Tutorials.html
## Ref: https://www.nature.com/articles/s41467-020-16904-3

##### Presetting ######
# rm(list = ls()) # Clean variable ##* Comment out if Run All
# memory.limit(150000)


#### Function: Run ROGUE ####
FUN_ROGUE <- function(seuratObject) {

  ## Load packages
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  if(!require("ROGUE")) devtools::install_github("PaulingLiu/ROGUE"); library(ROGUE)

  if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
  if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
  if(!require("dplyr")) install.packages("dplyr"); library(dplyr)

  # suppressMessages(library(ROGUE))
  # suppressMessages(library(ggplot2))
  # suppressMessages(library(tidyverse))

  expr <- seuratObject@assays[["RNA"]]@data %>% as.data.frame()
  meta <- seuratObject@meta.data

  ## Filtering out low-abundance genes and low-quality cells
  expr <- matr.filter(expr, min.cells = 10, min.genes = 10)

  ## Expression entropy model
  ent.res <- SE_fun(expr)

  # ## S-E plot
  # SEplot(ent.res)

  ## ROGUE calculation
  rogue.value <- CalculateRogue(ent.res, platform = "UMI")

  ## Calculate the ROGUE value of each putative cluster for each sample
  rogue.res <- rogue(expr, labels = meta$Cell_Type, samples = meta$Cell_Type, platform = "UMI", span = 0.6)
  # rogue.res <- rogue(expr, labels = meta$Cell_Type, samples = meta$Patient, platform = "UMI", span = 0.6)

  # ## Visualize ROGUE values on a boxplot
  # rogue.boxplot(rogue.res)

  Output.list <- list(ent.res = ent.res, rogue.value = rogue.value,
                      rogue.res = rogue.res)
  return(Output.list)
}


# ### Example usage:
# ## Load Demo
# load("D:/Dropbox/##_GitHub/###_VUMC/MagiCell/Demo/Seurat_PDAC_CombineObj_CTAnnot_BatchCor.RData")
#
# ## Store original cell type data
# seuratObject@meta.data$Cell_Type_Ori <- seuratObject@meta.data$Cell_Type
#
# ROGUE.list <- FUN_ROGUE(seuratObject)
#
# ## S-E plot
# SEplot(ROGUE.list[["ent.res"]])
#
# ## Visualize ROGUE values on a boxplot
# rogue.boxplot(ROGUE.list[["rogue.res"]])



#### Function: Recluster ####
FUN_Recluster <- function(seuratObject, rogue_res, ROGUE_Score = 0.9,
                          Num_PCA = 50, Cluster_Reso = 0.8) {
  # Iterate over cell types with scores < Set_ROGUE_Thr in rogue.res
  if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
  if(!require("dplyr")) install.packages("dplyr"); library(dplyr)

  # Ensure the 'CELL' column exists in the metadata, assuming it's equal to rownames if not already present
  if(!"CELL" %in% colnames(seuratObject@meta.data)) {
    seuratObject@meta.data$CELL <- rownames(seuratObject@meta.data)
  }

  for (cell_type in colnames(rogue_res)) {
    score <- rogue_res[cell_type, cell_type]


    if (!is.na(score) && score < ROGUE_Score) {
      # Subset the Seurat object for the specific cell type
      sub_seurat <- seuratObject[, seuratObject$Cell_Type == cell_type]


      # Perform clustering on the subset, adjust resolution as necessary
      sub_seurat <- FindNeighbors(sub_seurat, dims = 1:Num_PCA) %>%
        FindClusters(resolution = Cluster_Reso)


      # Create a dataframe with new labels
      new_cluster_labels <- data.frame(
        CELL = sub_seurat@meta.data$CELL,
        New_Cell_Type = paste(cell_type, sub_seurat$seurat_clusters, sep=" ")
      )


      # Join the new labels with the original Seurat object's metadata
      updated_metadata <- left_join(seuratObject@meta.data, new_cluster_labels, by = "CELL")


      # Update the Cell_Type in seuratObject with New_Cell_Type where not NA
      seuratObject@meta.data$Cell_Type <- ifelse(is.na(updated_metadata$New_Cell_Type),
                                                 seuratObject@meta.data$Cell_Type,
                                                 updated_metadata$New_Cell_Type)
    }
  }

  # The seuratObject now has updated Cell_Type with new subclusters for specific cell types

  return(seuratObject)
}

# ### Example usage:
# seuratObject <- FUN_Recluster(seuratObject, ROGUE.list[["rogue.res"]],Cluster_Reso = 0.6)
# DimPlot(seuratObject, reduction = "umap", group.by = "Cell_Type")


#********************************************************************************
#********************************************************************************
#### Function: Define Subclusters ####
# Define the function for ROGUE-based reclustering
Define_Subcluster <- function(seuratObject, ROGUE_Thr = 0.9, SubClust_Resolution = 0.8,
                              Cluster_Reso_Increment = 0.05, Num_PCA = 50) {

  # Load necessary packages
  if (!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
  if (!require("Seurat")) install.packages("Seurat"); library(Seurat)

  repeat {
    # Perform ROGUE analysis
    ROGUE.list <- FUN_ROGUE(seuratObject)

    # Check if any score is below ROGUE_Thr
    if (all(ROGUE.list[["rogue.res"]] >= ROGUE_Thr, na.rm = TRUE)) {
      break # Exit loop if all scores are >= ROGUE_Thr
    }

    # Re-cluster with the current resolution
    seuratObject <- FUN_Recluster(seuratObject, ROGUE.list[["rogue.res"]],
                                  Cluster_Reso = SubClust_Resolution,
                                  ROGUE_Score = ROGUE_Thr,
                                  Num_PCA = Num_PCA)

    # Decrease SubClust_Resolution for the next iteration
    SubClust_Resolution <- SubClust_Resolution + Cluster_Reso_Increment

    # Optionally, add a condition to prevent SubClust_Resolution from becoming too high
    if (SubClust_Resolution > 0.9) {
      warning("Resolution for FindClusters() has reached the threshold. Ending loop.")
      break
    }
  }

  return(seuratObject)
}

# ### Example usage:
# seuratObject <- Define_Subcluster(seuratObject, ROGUE_Thr = 0.9,
#                                   SubClust_Resolution = 0.6,
#                                   Num_PCA = 50)
#
#
# # After loop completion, seuratObject will have been re-clustered until the specified condition was met
#
# DimPlot(seuratObject, reduction = "umap", group.by = "Cell_Type")
# DimPlot(seuratObject, reduction = "umap", group.by = "seurat_clusters")


