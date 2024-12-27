##### Presetting #####
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Packages #####
if(!require('Seurat')) {install.packages('Seurat'); library(Seurat)}
if(!require('hdf5r')) {install.packages('hdf5r'); library(hdf5r)}
if(!require('Matrix')) {install.packages('Matrix'); library(Matrix)}
if(!require('tidyverse')) {install.packages('tidyverse'); library(tidyverse)}
if(!require('reshape2')) {install.packages('reshape2'); library(reshape2)}
if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}

##### Set Parameter #####
Set_GeneName <- c("EPGN", "EGFR") # 多個基因名稱
Set_DataType <- "GSE206391_ST"

# Set export parameters
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 10) # Generate unique time-based ID
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))
Set_note <- paste0(Name_FileID, "_", Set_DataType)
Name_Export <- paste0(Name_FileID)
Name_ExportFolder <- paste0("Export_", Set_note)
# Create export folder if it does not exist
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}


##### File Path #####
file_path <- "C:/Charlene/Dataset_Online/GSE206391/GSE206391_Preprocessed_data.h5"
h5_file <- H5File$new(file_path, mode = "r")

##### Extract Expression Data #####
counts_dataset <- h5_file[["layers/counts"]]
dims <- c(16685, 59319)  # Manual dimensions based on file inspection
expression_data <- counts_dataset[1:dims[1], 1:dims[2]]

# Extract gene and cell names
gene_names <- h5_file[["var/gene_name"]][]
cell_names <- h5_file[["obs/_index"]][]
rownames(expression_data) <- gene_names
colnames(expression_data) <- cell_names

# Convert to sparse matrix
sparse_matrix <- Matrix::Matrix(expression_data, sparse = TRUE)

##### Extract Metadata #####
obs_group <- h5_file[["obs"]]
metadata <- list()

for (dataset_name in names(obs_group)) {
  dataset <- obs_group[[dataset_name]]
  if (inherits(dataset, "H5D")) {
    metadata[[dataset_name]] <- dataset[]
  }
}

metadata <- as.data.frame(metadata)

# Check and process categories
if ("__categories" %in% names(obs_group)) {
  categories_group <- obs_group[["__categories"]]
  for (category_name in names(categories_group)) {
    if (category_name %in% colnames(metadata)) {
      category_mapping <- categories_group[[category_name]][]
      metadata[[category_name]] <- factor(metadata[[category_name]], levels = seq_along(category_mapping), labels = category_mapping)
    }
  }
}

##### Create Seurat Object #####
seurat_object <- CreateSeuratObject(counts = sparse_matrix, project = "GSE206391", meta.data = metadata)

##### Add Disease Labels #####
ad_patients <- c(2, 5, 8, 11, 15, 20, 34, 35, 36)
lp_patients <- c(3, 6, 9, 12, 14, 25, 26, 27, 28, 30, 37)
pso_patients <- c(1, 10, 13, 19, 22, 29, 31, 32, 33)

seurat_object$disease <- NA
seurat_object$disease[seurat_object$patient %in% ad_patients] <- "AD"
seurat_object$disease[seurat_object$patient %in% lp_patients] <- "LP"
seurat_object$disease[seurat_object$patient %in% pso_patients] <- "Pso"

##### Filter Valid Cells #####
valid_cells <- !is.na(seurat_object$patient)
seurat_object <- seurat_object[, valid_cells]

##### Add Group and GroupLN #####
seurat_object$group <- paste0(seurat_object$disease, "_", sub(".*_", "", seurat_object$replicates))
seurat_object$groupLN <- paste0(sub(".*_", "", seurat_object$replicates))

##### Preprocessing #####
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
seurat_object <- RunUMAP(seurat_object, dims = 1:10)

##### Visualize UMAP #####
DimPlot(seurat_object, reduction = "umap", group.by = "disease")
DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters")

##### Marker Analysis #####
# Set_GeneName <- c("EPGN", "EGFR", "KRT1", "KRT9", "KIT", "EP300", "NF1", "COL1A1")
VlnPlot(seurat_object, features = Set_GeneName[1:4], ncol = 2, group.by = "group")

##### Extract Data for Statistical Analysis #####
plot_data <- FetchData(seurat_object, vars = c(Set_GeneName, "group", "groupLN"))

##### Statistical Analysis #####
group_comparisons <- list(
  c("AD_NON LESIONAL", "AD_LESIONAL"),
  c("LP_NON LESIONAL", "LP_LESIONAL"),
  c("Pso_NON LESIONAL", "Pso_LESIONAL")
)

calculate_stats <- function(data, marker, groups, test_method = "t.test") {
  data_filtered <- data %>%
    filter(group %in% groups) %>%
    select(group, !!sym(marker)) %>%
    rename(Group = group, Expression = !!sym(marker))
  
  group_means <- data_filtered %>%
    group_by(Group) %>%
    summarize(mean_expression = mean(Expression), .groups = 'drop') %>%
    arrange(match(Group, groups))
  
  log2FC <- group_means$mean_expression[2] - group_means$mean_expression[1]
  FC <- 2^log2FC
  
  group1_values <- data_filtered %>% filter(Group == groups[1]) %>% pull(Expression)
  group2_values <- data_filtered %>% filter(Group == groups[2]) %>% pull(Expression)
  
  if (test_method == "t.test") {
    test_result <- t.test(group1_values, group2_values)
    p_value <- test_result$p.value
    statistic <- test_result$statistic
  } else if (test_method == "wilcox.test") {
    test_result <- wilcox.test(group1_values, group2_values)
    p_value <- test_result$p.value
    statistic <- test_result$statistic
  } else {
    stop("Unsupported test method: ", test_method)
  }
  
  return(data.frame(Marker = marker, Group1 = groups[1], Group2 = groups[2], log2FC = log2FC, FC = FC, p_value = p_value, statistic = statistic, test = test_method))
}

stats_results <- lapply(Set_GeneName, function(marker) {
  do.call(rbind, lapply(group_comparisons, function(groups) {
    calculate_stats(plot_data, marker, groups, test_method = "wilcox.test")
  }))
}) %>%
  do.call(rbind, .)

stats_results$p_adjusted <- p.adjust(stats_results$p_value, method = "BH")

rownames(stats_results) <- 1:nrow(stats_results)

# 將結果保存到文件
write.csv(stats_results, paste0(Name_ExportFolder, "/", Name_Export, "_stats.csv"), row.names = TRUE)

##### Boxplot Visualization #####
plot_data_long <- reshape2::melt(plot_data, id.vars = c("group", "groupLN"), variable.name = "Marker", value.name = "Expression")

ggplot(plot_data_long, aes(x = group, y = Expression, fill = groupLN)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.2, shape = 21, stroke = 0.5, color = "black") +
  facet_wrap(~ Marker, ncol = 2, scales = "free_y") +
  labs(x = "Group", y = "Expression Level", fill = "GroupLN") +
  theme_minimal(base_size = 15) +
  theme(
    strip.text = element_text(size = 24, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 22),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16)
  ) -> Plot_Box

Plot_Box

# 將圖形輸出到 PDF
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_BoxbPlot.pdf"), width = 12, height = 6)
print(Plot_Box)
dev.off()



##### Distribution Visualization #####
ggplot(plot_data_long, aes(x = Expression, fill = group)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ Marker, scales = "free", ncol = 2) +
  labs(title = "Expression Distribution by Marker", x = "Expression Level", y = "Density", fill = "Group") +
  theme_minimal(base_size = 15) +
  theme(
    strip.text = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )-> Plot_Distrib

Plot_Distrib

# 將圖形輸出到 PDF
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_DistribPlot.pdf"), width = 12, height = 6)
print(Plot_Distrib)
dev.off()


#### Export ####
## Export RData
save.image(paste0(Name_ExportFolder, "/", Name_Export, ".RData"))

