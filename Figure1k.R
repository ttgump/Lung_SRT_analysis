rm(list=ls())
library(Seurat)
library(stringr)
library(ggplot2)
library(CellChat)
library(patchwork)
library(pheatmap)
library(extrafont)
options('future.globals.maxSize'=1e16)

samples <- c("mock", "mild_day1", "mild_day7", "severe_day1", "severe_day7")
load("all_12_samples_Krt8.RData")

obj_list <- list()
for(i in 1:5) {
  obj <- readRDS(paste("../", samples[i], "_Seuratobject.rds", sep=""))
  
  load(paste("../", samples[i], "_RCTD_Plots/results.RData", sep=""))
  
  results_df <- results$results_df
  weight_doublet <- as.matrix(results$weights_doublet)
  colnames(weight_doublet) <- paste(colnames(weight_doublet), "score", sep="_")
  results_df$first_type <- as.character(results_df$first_type)
  results_df <- cbind(results_df, weight_doublet)
  results_df <- results_df[results_df$spot_class != "reject", ]
  results_df <- results_df[results_df$first_type != "VECs", ]
  
  results_df_selected <- results_df[results_df$first_type %in% c("AT1 cells", "AT2 cells", "Krt8 ADI", "Activated AT2 cells") & 
                                      results_df$first_type_score > 0.7, ]
  obj <- obj[, rownames(results_df_selected)]
  
  obj$cell_type <- as.character(results_df_selected$first_type)
  obj$sample <- samples[i]
  
  ### rename Transitional AT2 and ABMC
  colnames(obj) <- paste(samples[i], colnames(obj), sep="_")
  curr_krt8_obj <- krt8_obj[, krt8_obj$sample == samples[i]]
  obj$cell_type2 <- obj$cell_type
  obj$cell_type2[names(curr_krt8_obj$krt8_cluster)] <- as.character(curr_krt8_obj$krt8_cluster)
  
  obj <- obj[, obj$cell_type2 %in% c("AT1 cells", "AT2 cells", "Transitional AT2", "ABMC", "Activated AT2 cells")]
  
  obj_list[[i]] <- obj
}


merged_obj <- merge(obj_list[[1]], obj_list[2:5], add.cell.ids=samples)

merged_obj <- NormalizeData(merged_obj)
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])


data.input <- merged_obj[["RNA"]]$data
meta <- data.frame(labels = merged_obj$cell_type2, row.names = colnames(merged_obj))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 6) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type =  "truncatedMean", trim = 0.1)

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))

colors <- c("AT2"="#DEB887", "Activated AT2"="#6B8E23", "AT1"="lightblue",
            "Transitional AT2"="steelblue4", "ABMC"="#8B5A2B")


par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, color.use = colors, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, color.use = colors, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


ptm = Sys.time()
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "all")

ht_matrix <- ht@matrix
ht_matrix[is.na(ht_matrix)] <- 0


### plot Figure 1K
pheatmap(ht_matrix, color = viridis::viridis(100))
