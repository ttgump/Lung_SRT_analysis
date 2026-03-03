rm(list=ls())
library(Seurat)
library(SeuratExtend)
library(scater)
library(SingleCellExperiment)
library(stringr)
library(ggplot2)
library(dplyr)
library(umap)
library(RColorBrewer)
library(cowplot)
library(scCustomize)

### load AT cells (AT1, AT2, Activated AT2, Krt8 ADI) across 12 samples
samples <- c("mock", "mild_day1", "mild_day2", "mild_day3", "mild_day5", "mild_day7", "mild_day14",
             "severe_day1", "severe_day2", "severe_day3", "severe_day5", "severe_day7")

samples2 <- c("mock", "mild", "mild", "mild", "mild", "mild", "mild",
              "severe", "severe", "severe", "severe", "severe")

obj_list <- list()
for(i in 1:12) {
  obj <- readRDS(paste(samples[i], "_Seuratobject.rds", sep=""))
  
  load(paste(samples[i], "_RCTD_Results/results.RData", sep=""))
  
  results_df <- results$results_df
  weight_doublet <- as.matrix(results$weights_doublet)
  colnames(weight_doublet) <- paste(colnames(weight_doublet), "score", sep="_")
  results_df <- cbind(results_df, weight_doublet)
  results_df <- results_df[results_df$spot_class != "reject", ]
  results_df <- results_df[results_df$first_type != "VECs", ]
  
  results_df_selected <- results_df[results_df$first_type %in% c("AT1 cellss", "AT2 cells", "Activated AT2 cells", "Krt8 ADI") & 
                                      results_df$first_type_score > 0.7, ]
  obj <- obj[, rownames(results_df_selected)]
  
  obj$cell_type <- as.character(results_df_selected$first_type)
  obj$sample <- samples[i]
  obj$sample2 <- samples2[i]
  
  obj_list[[i]] <- obj
}

merged_obj <- merge(obj_list[[1]], obj_list[2:12], add.cell.ids=samples)
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)

merged_obj <- IntegrateLayers(object = merged_obj, method = JointPCAIntegration, orig.reduction = "pca", new.reduction = "jointpca",
                              verbose = TRUE)
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])
merged_obj <- FindNeighbors(merged_obj, dims = 1:20, reduction = "jointpca")
merged_obj <- RunUMAP(merged_obj, dims = 1:20, reduction = "jointpca")



emb_dat <- data.frame(merged_obj@reductions$umap@cell.embeddings, cell_type = merged_obj$cell_type, sample = merged_obj$sample, sample2 = merged_obj$sample2)
emb_dat$cell_type <- factor(emb_dat$cell_type, levels = c("AT2 cells", "Activated AT2 cells", "Krt8 ADI", "AT1 cellss"))
emb_dat$sample2 <- factor(emb_dat$sample2, levels = c("mock", "mild", "severe"))

p1 <- ggplot() +
  geom_point(data = emb_dat, mapping = aes(x = umap_1, y = umap_2, col = sample2), size=0.8) +
  labs(x="UMAP 1", y="UMAP 2", title = "Umap", col="Samples") +
  scale_color_manual(values = c("mock" = "#5F9EA0", "mild"="#6A5ACD", "severe"="#CD5C5C")) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

p2 <- ggplot() +
  geom_point(data = emb_dat, mapping = aes(x = umap_1, y = umap_2, col = cell_type), size=0.8) +
  labs(x="UMAP 1", y="UMAP 2", title = "Umap", col="cell_type") +
  scale_color_manual(values = c("AT2 cells"="#DEB887", "Activated AT2 cells"="#6B8E23", "AT1 cellss"="lightblue", "Krt8 ADI"="lightsalmon")) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# plot Figure 1E
plot_grid(p1, p2, scale=0.95)

### select Krt8 ADI cells and clustering to Transitional AT2 and ABMC
krt8_obj <- merged_obj[, merged_obj$cell_type == "Krt8 ADI"]
krt8_obj[["RNA"]] <- split(krt8_obj[["RNA"]], f = krt8_obj$sample)
krt8_obj <- NormalizeData(krt8_obj)
krt8_obj <- FindVariableFeatures(krt8_obj)
krt8_obj <- ScaleData(krt8_obj)
krt8_obj <- RunPCA(krt8_obj)

krt8_obj <- IntegrateLayers(object = krt8_obj, method = JointPCAIntegration, orig.reduction = "pca", k.weight = 50, new.reduction = "jointpca",
                            verbose = TRUE)
krt8_obj[["RNA"]] <- JoinLayers(krt8_obj[["RNA"]])
krt8_obj <- FindNeighbors(krt8_obj, dims = 1:20, reduction = "jointpca")
krt8_obj <- FindClusters(krt8_obj, resolution=0.4)
krt8_obj <- RunUMAP(krt8_obj, dims = 1:10, reduction = "jointpca")

DimPlot(krt8_obj, reduction = "umap", group.by = c("sample2", "seurat_clusters"), pt.size=1)

krt8_cluster <- rep("Krt8 ADI undetermined", ncol(krt8_obj))
krt8_cluster[krt8_obj$seurat_clusters %in% c(0, 2, 6, 9)] <- "Transitional AT2"
krt8_cluster[krt8_obj$seurat_clusters %in% c(1, 3, 4, 5, 8)] <- "ABMC"
names(krt8_cluster) <- colnames(krt8_obj)
krt8_cluster <- factor(krt8_cluster)
krt8_obj$krt8_cluster <- krt8_cluster


krt8_emb <- data.frame(krt8_obj@reductions$umap@cell.embeddings, cell_type = krt8_obj$krt8_cluster, sample = krt8_obj$sample, sample2 = krt8_obj$sample2)
krt8_emb$cell_type <- factor(krt8_emb$cell_type, levels = c("AT2 cells", "Activated AT2 cells", "AT1 cellss",  "Transitional AT2", "ABMC", "Krt8 ADI undetermined"))
krt8_emb$sample2 <- factor(krt8_emb$sample2, levels = c("mock", "mild", "severe"))

p1 <- ggplot() +
  geom_point(data = krt8_emb, mapping = aes(x = umap_1, y = umap_2, col = cell_type), size=0.8) +
  labs(x="UMAP 1", y="UMAP 2", title = "Umap", col="cell_type") +
  scale_color_manual(values = c("AT2 cells"="#DEB887", "Activated AT2 cells"="#6B8E23", "AT1 cellss"="lightblue", 
                                "Transitional AT2"="steelblue4", "ABMC"="#8B5A2B", "Krt8 ADI undetermined"="lightgrey")) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

p2 <- ggplot() +
  geom_point(data = krt8_emb, mapping = aes(x = umap_1, y = umap_2, col = sample2), size=0.8) +
  labs(x="UMAP 1", y="UMAP 2", title = "Umap", col="cell_type") +
  scale_color_manual(values = c("mock" = "#5F9EA0", "mild"="#6A5ACD", "severe"="#CD5C5C")) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# plot Figure 1F
plot_grid(p1, p2, scale=0.95, align="hv", nrow=1)

# plot Figure 1G
FeaturePlot_scCustom(krt8_obj, features = c("Krt5", "Krt17", "Itgb4", "Krt14"), keep.scale = NULL, colors_use = viridis_light_high, na_color = "lightgray", pt.size = 1.3)


### Differential analysis between Transitional AT2 and ABMC
krt8_obj@active.ident <- krt8_obj$krt8_cluster

krt8_obj_select <- krt8_obj[, krt8_obj$krt8_cluster != "Krt8 ADI undetermined"]

de_res <- FindAllMarkers(krt8_obj_select, logfc.threshold = 0, only.pos = T, min.pct = 0.1)


### rename cell types in all 12 samples
cell_type <- as.character(merged_obj$cell_type)
names(cell_type) <- colnames(merged_obj)
cell_type[colnames(krt8_obj)] <- as.character(krt8_obj$krt8_cluster)

merged_obj$cell_type_2 <- factor(cell_type, levels = rev(c("AT2 cells", "Activated AT2 cells", "AT1 cellss", "Transitional AT2", "ABMC", "Krt8 ADI undetermined")))
merged_obj <- merged_obj[, merged_obj$cell_type_2 != "ADI undetermined"]

# plot Figure 1H
DotPlot_scCustom(merged_obj, features = c("Sftpa1", "Sftpb", "Sftpc", "Abca3", "Lcn2", "Il33", "Ifi27l2a", "Lrg1",
                                          "Saa3", "Sprr1a", "Cdkn1a", "Ankrd1", "Ager", "Hopx", "Clic5", "Cav1",
                                          "Wfdc2", "Nupr1", "Clu", "Krt5", "Krt17"), group.by = "cell_type_2", flip_axes = F, x_lab_rotate = T, colors_use = rev(viridis_dark_high))


#### rename cell types in severe samples
merged_obj$cell_type2 <- merged_obj$cell_type
merged_obj$cell_type2[names(krt8_obj$krt8_cluster)] <- as.character(krt8_obj$krt8_cluster)

severe_obj <- merged_obj[, merged_obj$sample2 %in% c("mock", "severe")]
severe_obj[["RNA"]] <- split(severe_obj[["RNA"]], f = severe_obj$sample)
severe_obj <- NormalizeData(severe_obj)
severe_obj <- FindVariableFeatures(severe_obj)
severe_obj <- ScaleData(severe_obj)
severe_obj <- RunPCA(severe_obj)

severe_obj <- IntegrateLayers(object = severe_obj, method = JointPCAIntegration, orig.reduction = "pca", new.reduction = "jointpca",
                            verbose = TRUE)
severe_obj[["RNA"]] <- JoinLayers(severe_obj[["RNA"]])
severe_obj <- FindNeighbors(severe_obj, dims = 1:20, reduction = "jointpca")
severe_obj <- RunUMAP(severe_obj, dims = 1:20, reduction = "jointpca")


severe_emb <- data.frame(severe_obj@reductions$umap@cell.embeddings, cell_type = severe_obj$cell_type2, sample = severe_obj$sample, sample2 = severe_obj$sample2)
severe_emb$cell_type <- factor(severe_emb$cell_type, levels = c("AT2 cells", "Activated AT2 cells", "AT1 cellss", "Transitional AT2", "ABMC", "Krt8 ADI undetermined"))

### Plot Figure 1I
ggplot() +
  geom_point(data = severe_emb, mapping = aes(x = umap_1, y = umap_2, col = cell_type), size=0.8) +
  labs(x="UMAP 1", y="UMAP 2", title = "Umap", col="cell_type") +
  scale_color_manual(values = c("AT2 cells"="#DEB887", "Activated AT2 cells"="#6B8E23", "AT1 cellss"="lightblue", "Krt8 ADI"="lightsalmon",
                                "Transitional AT2"="steelblue4", "ABMC"="#8B5A2B", "Krt8 ADI undetermined"="lightgrey")) +
  facet_wrap(~sample, nrow = 1, scales = "free") + theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))


severe_obj2 <- severe_obj[severe_obj$sample %in% c("mock", "zhong_day1", "zhong_day2", "zhong_day7"), ]
severe_obj_agg2 <- aggregate(severe_obj2$count, by = list(severe_obj2$celltype, severe_obj2$sample), sum)
colnames(severe_obj_agg2) <- c("celltype", "sample", "Freq")
severe_obj_agg2$celltype <- factor(severe_obj_agg2$celltype, levels = c("AT1", "AT2", "Activated AT2", "Transitional AT2", "ABMC"))
severe_obj_agg2$sample <- factor(severe_obj_agg2$sample, levels = c("mock", "zhong_day1", "zhong_day2", "zhong_day7"))

ggplot(severe_obj_agg2) +
  aes(x = sample, fill = celltype, weight = Freq, by = sample) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("AT2"="#DEB887", "Activated AT2"="#6B8E23", "AT1"="lightblue",
                               "Transitional AT2"="steelblue4", "ABMC"="#8B5A2B", "Krt8 ADI undetermined"="lightgrey")) +
  theme_classic()


#### rename cell types in mild samples
mild_obj <- merged_obj[, merged_obj$sample2 %in% c("mock", "mild")]
mild_obj[["RNA"]] <- split(mild_obj[["RNA"]], f = mild_obj$sample)
mild_obj <- NormalizeData(mild_obj)
mild_obj <- FindVariableFeatures(mild_obj)
mild_obj <- ScaleData(mild_obj)
mild_obj <- RunPCA(mild_obj)

mild_obj <- IntegrateLayers(object = mild_obj, method = JointPCAIntegration, orig.reduction = "pca", new.reduction = "jointpca",
                             verbose = TRUE)
mild_obj[["RNA"]] <- JoinLayers(mild_obj[["RNA"]])
mild_obj <- FindNeighbors(mild_obj, dims = 1:20, reduction = "jointpca")
mild_obj <- RunUMAP(mild_obj, dims = 1:20, reduction = "jointpca")


mild_emb <- data.frame(mild_obj@reductions$umap@cell.embeddings, cell_type = mild_obj$cell_type2, sample = mild_obj$sample, sample2 = mild_obj$sample2)
mild_emb$cell_type <- factor(mild_emb$cell_type, levels = c("AT2 cells", "Activated AT2 cells", "AT1 cellss", "Transitional AT2", "ABMC", "Krt8 ADI undetermined"))


ggplot() +
  geom_point(data = mild_emb, mapping = aes(x = -umap_1, y = umap_2, col = cell_type), size=0.8) +
  labs(x="UMAP 1", y="UMAP 2", title = "Umap", col="cell_type") +
  scale_color_manual(values = c("AT2 cells"="#DEB887", "Activated AT2 cells"="#6B8E23", "AT1 cellss"="lightblue", "Krt8 ADI"="lightsalmon",
                                "Transitional AT2"="steelblue4", "ABMC"="#8B5A2B", "Krt8 ADI undetermined"="lightgrey")) +
  facet_wrap(~sample, nrow = 1, scales = "free") + theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))
