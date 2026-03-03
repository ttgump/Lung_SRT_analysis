rm(list=ls())
library(Seurat)
library(SeuratExtend)
library(slingshot)
library(stringr)
library(ggplot2)
library(anndataR)
library(dplyr)
library(umap)
library(RColorBrewer)
library(mclust)
library(umap)
library(forcats)
library(cowplot)
library(scCustomize)


load("mock_severes_ATs.RData")
severe_emb_dat <- data.frame(severe_obj@reductions$umap@cell.embeddings, cell_type = severe_obj$cell_type2, sample = severe_obj$sample)
severe_emb_dat <- severe_emb_dat[severe_emb_dat$cell_type != "Krt8 ADI undetermined", ]

ggplot() +
  geom_point(data = severe_emb_dat, mapping = aes(x = umap_1, y = umap_2, col = cell_type), size=0.5) +
  labs(x="UMAP 1", y="UMAP 2", title = "Slingshot trajectory", col="Cluster") +
  scale_color_manual(values = c("AT2 cells"="#DEB887", "Activated AT2 cells"="#6B8E23", "AT1 cellss"="lightblue",
                                "Transitional AT2"="steelblue4", "AMBC"="#8B5A2B", "Krt8 ADI undetermined"="lightgrey")) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))


# mclust of the umap embedding
mclust <- Mclust(severe_emb_dat[, 1:2], G=6)
severe_emb_dat$G <-factor(mclust$classification)

ggplot() +
  geom_point(data = severe_emb_dat, mapping = aes(x = umap_1, y = umap_2, col = G), size=0.5) +
  labs(x="UMAP 1", y="UMAP 2", title = "Slingshot trajectory", col="Cluster") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))


### slingshot analysis
set.seed(1)
lineages <- getLineages(data = severe_emb_dat[, 1:2], clusterLabels = severe_emb_dat$G, start.clus = 4, end.clus = c(2, 6))
lineages
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves

severe_emb_dat$cell_type <- factor(severe_emb_dat$cell_type, levels = c("AT2 cells", "Activated AT2 cells", "AT1 cellss", "Transitional AT2", "ABMC"))


p1 <- ggplot() +
  geom_point(data = severe_emb_dat, mapping = aes(x = umap_1, y = umap_2, col = cell_type), size=0.5) +
  scale_color_manual(values = c("AT2 cells"="#DEB887", "Activated AT2 cells"="#6B8E23", "AT1 cellss"="lightblue",
                                "Transitional AT2"="steelblue4", "AMBC"="#8B5A2B", "Krt8 ADI undetermined"="lightgrey")) +
  labs(x="UMAP 1", y="UMAP 2", title = "Slingshot trajectory", col="Cluster") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size = 3)))



severe_emb_dat$pseudotime <- apply(curves@assays@data$pseudotime, 1, function(z) {return(as.numeric(mean(z, na.rm=T)))})
p2 <- ggplot() +
  geom_point(data = severe_emb_dat, mapping = aes(x = umap_1, y = umap_2, col = pseudotime), size=0.5) +
  scale_colour_gradientn(colors = viridis::viridis(n=50)) +
  labs(x="UMAP 1", y="UMAP 2", title = "Slingshot trajectory", col="pseudotime") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5))

# Adding the curves
for (i in seq_along(slingCurves(curves))) {
  curve_i <- slingCurves(curves)[[i]]
  curve_i <- curve_i$s[curve_i$ord, ]
  colnames(curve_i) <- c("UMAP1", "UMAP2")
  p2 <- p2 + geom_path(data = as.data.frame(curve_i), 
                       aes(x = -UMAP1, y = UMAP2), linewidth=1)
}

### plot Figure 1J
plot_grid(p1, p2, align = "hv", scale = 0.95, nrow = 1)


ggplot(severe_emb_dat, aes(x = pseudotime , fill = cell_type)) +
  facet_wrap(~sample, ncol = 1, scales = "fixed", strip.position = "left") +
  geom_density(alpha = 0.75) + 
  scale_fill_manual(values = c("AT2 cells"="#DEB887", "Activated AT2 cells"="#6B8E23", "AT1 cellss"="lightblue",
                                "Transitional AT2"="steelblue4", "AMBC"="#8B5A2B", "Krt8 ADI undetermined"="lightgrey")) +
  theme_classic()

