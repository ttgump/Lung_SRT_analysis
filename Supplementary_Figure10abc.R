rm(list=ls())
library(Seurat)
library(scCustomize)
library(stringr)
library(ggplot2)
library(CellChat)
library(patchwork)
options('future.globals.maxSize'=1e16)


samples <- c("severe_day1", "severe_day2", "severe_day3", "severe_day5", "severe_day7")


### load data
obj_list <- list()

load("../all_12_samples_Krt8.RData")

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
  
  results_df$first_type[results_df$first_type %in% c("Recruited macrophages", "M2 macrophages")] <- "MDM"
  results_df$first_type[results_df$first_type %in% c("AM (PBS)", "Cd163- Cd11c+ IMs", "Cd163+ Cd11c+ IMs")] <- "AM + IM"
  
  results_df_selected <- results_df[results_df$first_type %in% c("AT1 cellss", "AT2 cells", "Krt8 ADI", "AM + IM", "B-lymphocytes", "T-lymphocytes", "MDM") & 
                                      results_df$first_type_score > 0.7, ]
  obj <- obj[, rownames(results_df_selected)]
  
  obj$cell_type <- as.character(results_df_selected$first_type)
  obj$sample <- samples[i]
  
  colnames(obj) <- paste(samples[i], colnames(obj), sep="_")
  curr_krt8_obj <- krt8_obj[, krt8_obj$sample == samples[i]]
  obj$cell_type2 <- obj$cell_type
  obj$cell_type2[names(curr_krt8_obj$krt8_cluster)] <- as.character(curr_krt8_obj$krt8_cluster)
  
  obj <- obj[, obj$cell_type2 %in% c("AT1 cellss", "AT2 cells", "Krt8 ADI c2", "AM + IM", "B-lymphocytes", "T-lymphocytes", "MDM")]
  
  obj_list[[i]] <- obj
}

### merge severe samples
merged_obj <- merge(obj_list[[1]], obj_list[2:5], add.cell.ids=samples)

merged_obj <- NormalizeData(merged_obj)
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])


data.input <- merged_obj[["RNA"]]$data
meta <- data.frame(labels = merged_obj$cell_type2, row.names = colnames(merged_obj))
meta$labels[meta$labels == "AT1 cellss"] <- "AT1"
meta$labels[meta$labels == "AT2 cells"] <- "AT2"
meta$labels[meta$labels == "Krt8 ADI c2"] <- "ABMC"
meta$labels[meta$labels == "AM + IM"] <- "AM"

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


pdf("5severe_cellchat.pdf", width=10, height=5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

save(cellchat, file="5severe_cellChat_objs.RData")
load("5severe_cellChat_objs.RData")


mat <- cellchat@net$weight

### plot communications for each cell type
for(i in 1:7) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  pdf(paste(rownames(mat)[i], "_5severe_cellchat.pdf", sep=""), width=5, height=5)
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}


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
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 20)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 20)
pdf("5severe_cellchat_signalingRole_heatmap.pdf", width=16, height=24)
ht1 + ht2
dev.off()


### plot outgoing and incoming patterns
library(NMF)
library(ggalluvial)

selectK(cellchat, pattern = "outgoing")

nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
pdf("5severe_cellchat_outgoing_river_plot.pdf", width=10, height=10)
netAnalysis_river(cellchat, pattern = "outgoing", cutoff = 0.4)
dev.off()
# dot plot
pdf("5severe_cellchat_outgoing_dot_plot.pdf", width=16, height=8)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()





selectK(cellchat, pattern = "incoming")

nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
pdf("5severe_cellchat_incoming_river_plot.pdf", width=10, height=10)
netAnalysis_river(cellchat, pattern = "incoming", cutoff = 0.4)
dev.off()
# dot plot
pdf("5severe_cellchat_incoming_dot_plot.pdf", width=16, height=8)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()
