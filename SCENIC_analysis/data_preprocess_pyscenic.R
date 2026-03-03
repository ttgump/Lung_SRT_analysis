library(SCENIC)
library(AUCell)
library(RcisTarget)
library(Seurat)

# Set working directory
setwd("/home/nas2/biod/luyuejing/R/TF_target")

# --- 1. Preprocessing Mock Dataset ---
obj_mock <- readRDS("./data/mock_Seuratobject.rds")
load("./data/mock_myRCTD_results.RData")

# Integrate RCTD doublet scores
results_df_mock <- results$results_df
weight_doublet_mock <- as.matrix(results$weights_doublet)
colnames(weight_doublet_mock) <- paste(colnames(weight_doublet_mock), "score", sep = "_")
results_df_mock <- cbind(results_df_mock, weight_doublet_mock)

# Filter low quality spots and standardize cell types (Merge Macrophages -> MDM)
results_df_mock <- results_df_mock[results_df_mock$spot_class != "reject", ]
results_df_replace_mock = results_df_mock
results_df_replace_mock$first_type <- as.character(results_df_replace_mock$first_type)
results_df_replace_mock$first_type[results_df_replace_mock$first_type %in% c("Recruited macrophages", "M2 macrophages")] <- "MDM"

# Select high-confidence cells for specific cell types
results_df_selected_mock <- results_df_replace_mock[results_df_replace_mock$first_type %in% c("AT1 cellss", "AT2 cells","AM (PBS)","MDM","B-lymphocytes","T-lymphocytes")& results_df_replace_mock$first_type_score > 0.6, ]
results_df_replace_mock$first_type <- as.factor(results_df_replace_mock$first_type)

# Update Seurat object
obj_mock <- obj_mock[, rownames(results_df_selected_mock)]
obj_mock[["cell_type"]] = as.factor(results_df_selected_mock$first_type)
results_df_selected_mock$first_type = as.factor(results_df_selected_mock$first_type)
celltypes_mock <- levels(results_df_selected_mock$first_type)

# --- 2. Preprocessing Qing_day1 Dataset ---
obj_qing_day1 <- readRDS("./data/qing_day1_Seuratobject.rds")
load("./data/qing_day1_myRCTD_results.RData")

# Integrate results and filter
results_df_qing_day1 <- results$results_df
weight_doublet_qing_day1 <- as.matrix(results$weights_doublet)
colnames(weight_doublet_qing_day1) <- paste(colnames(weight_doublet_qing_day1), "score", sep = "_")
results_df_qing_day1 <- cbind(results_df_qing_day1, weight_doublet_qing_day1)
results_df_qing_day1 <- results_df_qing_day1[results_df_qing_day1$spot_class != "reject", ]

# Standardize cell types
results_df_replace_qing_day1 = results_df_qing_day1
results_df_replace_qing_day1$first_type <- as.character(results_df_replace_qing_day1$first_type)
results_df_replace_qing_day1$first_type[results_df_replace_qing_day1$first_type %in% c("Recruited macrophages", "M2 macrophages")] <- "MDM"

# Select high-confidence cells
results_df_selected_qing_day1 <- results_df_replace_qing_day1[results_df_replace_qing_day1$first_type %in% c("AT1 cellss", "AT2 cells","AM (PBS)","MDM","B-lymphocytes","T-lymphocytes")& results_df_replace_qing_day1$first_type_score > 0.6, ]
results_df_replace_qing_day1$first_type <- as.factor(results_df_replace_qing_day1$first_type)

# Update Seurat object
obj_qing_day1 <- obj_qing_day1[, rownames(results_df_selected_qing_day1)]
obj_qing_day1[["cell_type"]] = as.factor(results_df_selected_qing_day1$first_type)
results_df_selected_qing_day1$first_type = as.factor(results_df_selected_qing_day1$first_type)
celltypes_qing_day1 <- levels(results_df_selected_qing_day1$first_type)

# --- 3. Preprocessing Zhong_day1 Dataset ---
obj_zhong_day1 <- readRDS("./data/zhong_day1_Seuratobject.rds")
load("./data/zhong_day1_myRCTD_results.RData")

# Integrate results and filter
results_df_zhong_day1 <- results$results_df
weight_doublet_zhong_day1 <- as.matrix(results$weights_doublet)
colnames(weight_doublet_zhong_day1) <- paste(colnames(weight_doublet_zhong_day1), "score", sep = "_")
results_df_zhong_day1 <- cbind(results_df_zhong_day1, weight_doublet_zhong_day1)
results_df_zhong_day1 <- results_df_zhong_day1[results_df_zhong_day1$spot_class != "reject", ]

# Standardize cell types
results_df_replace_zhong_day1 = results_df_zhong_day1
results_df_replace_zhong_day1$first_type <- as.character(results_df_replace_zhong_day1$first_type)
results_df_replace_zhong_day1$first_type[results_df_replace_zhong_day1$first_type %in% c("Recruited macrophages", "M2 macrophages")] <- "MDM"

# Select high-confidence cells
results_df_selected_zhong_day1 <- results_df_replace_zhong_day1[results_df_replace_zhong_day1$first_type %in% c("AT1 cellss", "AT2 cells","AM (PBS)","MDM","B-lymphocytes","T-lymphocytes")& results_df_replace_zhong_day1$first_type_score > 0.6, ]
results_df_replace_zhong_day1$first_type <- as.factor(results_df_replace_zhong_day1$first_type)

# Update Seurat object
obj_zhong_day1 <- obj_zhong_day1[, rownames(results_df_selected_zhong_day1)]
obj_zhong_day1[["cell_type"]] = as.factor(results_df_selected_zhong_day1$first_type)
results_df_selected_zhong_day1$first_type = as.factor(results_df_selected_zhong_day1$first_type)
celltypes_zhong_day1 <- levels(results_df_selected_zhong_day1$first_type)

# --- 4. Feature Selection & Merging ---

# Identify common Highly Variable Genes (HVGs) across datasets
data1_ <- FindVariableFeatures(obj_mock, nfeatures = 6000)
data2_ <- FindVariableFeatures(obj_qing_day1, nfeatures = 6000)
data3_ <- FindVariableFeatures(obj_zhong_day1, nfeatures = 6000)

hv_genes1 <- VariableFeatures(data1_)
hv_genes2 <- VariableFeatures(data2_)
hv_genes3 <- VariableFeatures(data3_)
# hv_genes1 <- rownames(data1) # Option to use all genes
# ...

common_hv_genes <- intersect(intersect(hv_genes1, hv_genes2), hv_genes3)

# Add specific genes of interest
additional_genes <- c("Nfkb1", "Hif1a", "Ppargc1a", "Ppara", "Pparg", "Ppard")
final_gene_set <- unique(c(common_hv_genes, additional_genes))

# Subset datasets to selected features
data1_subset <- obj_mock[final_gene_set]
data2_subset <- obj_qing_day1[final_gene_set]
data3_subset <- obj_zhong_day1[final_gene_set]

# Extract count matrices
matrix1 = GetAssayData(object = data1_subset, layer = "counts")
matrix2 = GetAssayData(object = data2_subset, layer = "counts")
matrix3 = GetAssayData(object = data3_subset, layer = "counts")


# Merge matrices and metadata
combined_matrix <- cbind(matrix1, matrix2, matrix3)

cluster1 = obj_mock[["cell_type"]] 
cluster2 = obj_qing_day1[["cell_type"]] 
cluster3 = obj_zhong_day1[["cell_type"]] 
cell_type_all = rbind(cluster1, cluster2, cluster3)

# --- 5. Export Data for External Tools ---
# Save metadata and expression matrix (Matrix Market format)
# Save cell metadata

write.csv(cell_type_all, "./filtered_feature_bc_matrix/cell_type_all_0.6_6000.csv", row.names = TRUE, col.names = TRUE, sep = ",", quote = TRUE)

library(Matrix)
# Save sparse matrix (Matrix Market format)
writeMM(combined_matrix, "./filtered_feature_bc_matrix/matrix_all_0.6_6000.mtx")

# Save gene names
write.table(
  data.frame(gene = rownames(combined_matrix)),
  file = "./filtered_feature_bc_matrix/genes_all_0.6_6000.tsv",
  row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
)

# Save barcodes
write.table(
  data.frame(cell = colnames(combined_matrix)),
  file = "./filtered_feature_bc_matrix/barcodes_all_0.6_6000.tsv",
  row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
)


