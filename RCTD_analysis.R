rm(list=ls())
library(Seurat)
library(SeuratExtend)
library(spacexr)
library(stringr)

### load high resolution reference data
load("mouse_references/EpiHiRes_seurat.RData")
EpiHiRes.obj <- UpdateSeuratObject(subset)
EpiHiRes.meta <- read.table("mouse_references/GSE141259_HighResolution_cellinfo.csv", sep="\t",
                            header=T, row.names=1)

#### rename cell types
EpiHiRes.meta <- EpiHiRes.meta[EpiHiRes.meta$cell_type %in% c("AT2", "Goblet", "Club", "AT1", "Ciliated", "Ciliated activated",
                                                              "Mki67+ Proliferation", "AT2 activated", "Krt8+ ADI", "Club to ciliated"), ]
EpiHiRes.meta$cell_type[EpiHiRes.meta$cell_type=="AT2"] <- "AT2 cells"
EpiHiRes.meta$cell_type[EpiHiRes.meta$cell_type=="Goblet"] <- "Goblet cells"
EpiHiRes.meta$cell_type[EpiHiRes.meta$cell_type=="Club"] <- "Club cells"
EpiHiRes.meta$cell_type[EpiHiRes.meta$cell_type=="AT1"] <- "AT1 cellss"
EpiHiRes.meta$cell_type[EpiHiRes.meta$cell_type=="Ciliated"] <- "Ciliated cells"
EpiHiRes.meta$cell_type[EpiHiRes.meta$cell_type=="Ciliated activated"] <- "Ciliated cell subset"
EpiHiRes.meta$cell_type[EpiHiRes.meta$cell_type=="Mki67+ Proliferation"] <- "Mki67+ proliferating cells"
EpiHiRes.meta$cell_type[EpiHiRes.meta$cell_type=="AT2 activated"] <- "Activated AT2 cells"
EpiHiRes.meta$cell_type[EpiHiRes.meta$cell_type=="Krt8+ ADI"] <- "Krt8 ADI"
common.cell1 <- intersect(colnames(EpiHiRes.obj), rownames(EpiHiRes.meta))
EpiHiRes.obj <- EpiHiRes.obj[, common.cell1]
EpiHiRes.meta <- EpiHiRes.meta[common.cell1, ]
EpiHiRes.obj$cell_type <- EpiHiRes.meta$cell_type
 
 ### load whole lung reference data
load("mouse_references/Whole_lung_seurat.RData")
WholeLung.obj <- UpdateSeuratObject(seu.ica)
WholeLung.meta <- read.csv("mouse_references//GSE141259_WholeLung_cellinfo.csv", stringsAsFactors=F, row.names=1)
WholeLung.meta <- WholeLung.meta[! WholeLung.meta$cell.type %in% c("AM (Bleo)", "CECs", 
                                                                   "Low quality cells", "LECs", "Fn1+ macrophages", "T cell subset",
                                                                   "Themis+ T-lymphocytes"), ]
common.cell2 <- intersect(colnames(WholeLung.obj), rownames(WholeLung.meta))
WholeLung.obj <- WholeLung.obj[, common.cell2]
WholeLung.meta <- WholeLung.meta[common.cell2, ]
WholeLung.obj$cell_type <- WholeLung.meta$cell.type

combined.obj <- merge(EpiHiRes.obj, WholeLung.obj, collapse=T, project="combined_obj")
combined.celltypes <- combined.obj$cell_type
combined.cell.idx <- which(combined.celltypes %in% names(table(combined.celltypes))[table(combined.celltypes) > 200] &
                             !is.na(combined.celltypes))
combined.obj <- combined.obj[, combined.cell.idx]


combined.obj.counts <- combined.obj@assays$RNA$counts
cell_types <- combined.obj$cell_type; names(cell_types) <- colnames(combined.obj)
cell_types[cell_types=="Cd163-/Cd11c+ IMs"] <- "Cd163- Cd11c+ IMs"
cell_types[cell_types=="Cd163+/Cd11c- IMs"] <- "Cd163+ Cd11c+ IMs"
cell_types[cell_types=="Mki67+/Top2a+ proliferating cells"] <- "Mki67+ Top2a+ proliferating cells"
cell_types <- as.factor(cell_types)
reference <- Reference(combined.obj.counts, cell_types)
print(dim(reference@counts))
print(table(reference@cell_types))

samples <- c("mock", "qing_day1", "qing_day2", "qing_day3", "qing_day5", "qing_day7", "qing_day14",
             "zhong_day1", "zhong_day2", "zhong_day3", "zhong_day5", "zhong_day7")


for(curr_sample in samples) {
  rna.data <- readRDS(paste(curr_sample, "_Seuratobject.rds", sep=""))
  barcode <- colnames(rna.data)
  barcode <- substr(barcode, 1, 19)
  barcode <- str_split(barcode, "_", simplify = T)
  spot_x <- as.integer(barcode[, 3])
  spot_y <- as.integer(barcode[, 4])
  spot_pos <- data.frame(x=spot_x, y=spot_y, row.names=colnames(rna.data))

  puck <- SpatialRNA(spot_pos, rna.data@assays$RNA$counts)
  print(dim(puck@counts)) # observe Digital Gene Expression matrix
  print(head(puck@coords)) # start of coordinate data.frame

  myRCTD <- create.RCTD(puck, reference, max_cores = 20)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')


  results <- myRCTD@results
  # normalize the cell type proportions to sum to 1.
  norm_weights = normalize_weights(results$weights) 
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- myRCTD@spatialRNA
  resultsdir <- paste(curr_sample, "_RCTD_Results", sep="") ## you may change this to a more accessible directory on your computer.
  dir.create(resultsdir)

  save(results, file=paste(curr_sample, "_RCTD_Results/results.RData", sep=""))

  plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)

  plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 

  plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                     results$results_df) 

  plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

  plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir) 

  doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",] 

  plot_doublets(spatialRNA, doublets, resultsdir, cell_type_names) 

  plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names) 

  doub_occur <- table(doublets$second_type, doublets$first_type) 

  plot_doub_occur_stack(doub_occur, resultsdir, cell_type_names) 
}
