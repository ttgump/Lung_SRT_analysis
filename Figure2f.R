rm(list=ls())
library(Seurat)
library(scCustomize)
library(stringr)
library(ggplot2)

samples <- c("mock", "mild_day1", "severe_day1")

samples2 <- c("mock", "mild", "severe")

obj_list <- list()
for(i in 1:3) {
  obj <- readRDS(paste("../", samples[i], "_Seuratobject.rds", sep=""))
  
  load(paste("../", samples[i], "_RCTD_Results/results.RData", sep=""))
  
  results_df <- results$results_df
  weight_doublet <- as.matrix(results$weights_doublet)
  colnames(weight_doublet) <- paste(colnames(weight_doublet), "score", sep="_")
  results_df$first_type <- as.character(results_df$first_type)
  results_df <- cbind(results_df, weight_doublet)
  results_df <- results_df[results_df$spot_class != "reject", ]
  results_df <- results_df[results_df$first_type != "VECs", ]
  results_df$first_type[results_df$first_type %in% c("Recruited macrophages", "M2 macrophages")] <- "MDM"
  
  results_df_selected <- results_df[results_df$first_type %in% c("AT1 cells", "AT2 cells", "AM (PBS)", "B-lymphocytes", "T-lymphocytes", "MDM"), ]
  obj <- obj[, rownames(results_df_selected)]
  
  obj$cell_type <- as.character(results_df_selected$first_type)
  obj$sample <- samples[i]
  obj$sample2 <- samples2[i]
  
  obj_list[[i]] <- obj
}

merged_obj <- merge(obj_list[[1]], obj_list[2:3], add.cell.ids=samples)

merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])

celltypes <- c("AT1 cells", "AT2 cells", "AM (PBS)", "B-lymphocytes", "T-lymphocytes", "MDM")
celltypes2 <- c("AT1", "AT2", "AM", "BLM", "TLM", "MDM")
genename <- c("Nfkb1", "Hif1a", "Ppargc1a", "Ppara", "Pparg", "Ppard")


### load regulon list across samples
gene_list <- list()
for(i in 1:6) {
  regulon <- read.csv(paste("results_6000/results_TF_target_", celltypes2[i], "/regulon_mat_all_regulons.csv", sep = ""), row.names = 1)
  
  curr_gene_list <- list()
  for(j in 1:6) {
    regulon_select <- regulon[rownames(regulon) %in% c(paste(genename[j], "(+)", sep=""), paste(genename[j], "(-)", sep="")), ]
    
    genename_2 <- apply(regulon_select, 1, function(z) colnames(regulon_select)[which(z==1)])
    genename_2 <- unlist(genename_2)
    genename_2 <- as.character(genename_2)
    genename_2 <- gsub("\\.", "-", genename_2)
    genename_2 <- genename_2[!duplicated(genename_2)]
    
    curr_gene_list[[j]] <- as.character(genename_2)
  }
  names(curr_gene_list) <- genename
  
  gene_list[[i]] <- curr_gene_list
}
names(gene_list) <- celltypes


### aggregate each gene's regulon and remove duplicates
gene_list_2 <- list()
for(i in 1:6) {
  curr_gene_list <- list()
  for(j in 1:6) {
    curr_gene_list[[j]] <- gene_list[[celltypes[j]]][[genename[i]]]
  }
  curr_gene_list <- unlist(curr_gene_list)
  curr_gene_list <- curr_gene_list[!duplicated(curr_gene_list)]
  gene_list_2[[i]] <- curr_gene_list
}
names(gene_list_2) <- genename



### calculate module scores
gene_list_2 <- gene_list_2[-3]
merged_obj <- AddModuleScore(merged_obj, features = gene_list_2, name=names(gene_list_2))
merged_obj$celltype_stage <- paste(merged_obj$cell_type, merged_obj$sample, sep="-")
merged_obj$celltype_stage <- factor(merged_obj$celltype_stage, levels = c("AT1 cells-severe_day1", "AT2 cells-severe_day1", "AM (PBS)-severe_day1", "B-lymphocytes-severe_day1", "T-lymphocytes-severe_day1", "MDM-severe_day1",
                                                                          "AT1 cells-mild_day1", "AT2 cells-mild_day1", "AM (PBS)-mild_day1", "B-lymphocytes-mild_day1", "T-lymphocytes-mild_day1", "MDM-mild_day1",
                                                                          "AT1 cells-mock", "AT2 cells-mock", "AM (PBS)-mock", "B-lymphocytes-mock", "T-lymphocytes-mock", "MDM-mock"))

pdf("Nfkb1_target_genes_plot.pdf", width=13, height=13.5, onefile=F)
RidgePlot(merged_obj,features=c('Nfkb11'), group.by='celltype_stage', fill.by="sample", cols = c("AT1 cells-severe_day1"="#CD5C5C", "AT2 cells-severe_day1"="#CD5C5C", "AM (PBS)-severe_day1"="#CD5C5C", 
                                                                                                 "B-lymphocytes-severe_day1"="#CD5C5C", "T-lymphocytes-severe_day1"="#CD5C5C", "MDM-severe_day1"="#CD5C5C",
                                                                                                 "AT1 cells-mild_day1"="lightblue", "AT2 cells-mild_day1"="lightblue", "AM (PBS)-mild_day1"="lightblue", 
                                                                                                 "B-lymphocytes-mild_day1"="lightblue", "T-lymphocytes-mild_day1"="lightblue", "MDM-mild_day1"="lightblue",
                                                                                                 "AT1 cells-mock"="lightgrey", "AT2 cells-mock"="lightgrey", "AM (PBS)-mock"="lightgrey", 
                                                                                                 "B-lymphocytes-mock"="lightgrey", "T-lymphocytes-mock"="lightgrey", "MDM-mock"="lightgrey"))
dev.off()

pdf("Hif1a_target_genes_plot.pdf", width=13, height=13.5, onefile=F)
RidgePlot(merged_obj,features=c('Hif1a2'), group.by='celltype_stage', fill.by="sample", cols = c("AT1 cells-severe_day1"="#CD5C5C", "AT2 cells-severe_day1"="#CD5C5C", "AM (PBS)-severe_day1"="#CD5C5C", 
                                                                                                 "B-lymphocytes-severe_day1"="#CD5C5C", "T-lymphocytes-severe_day1"="#CD5C5C", "MDM-severe_day1"="#CD5C5C",
                                                                                                 "AT1 cells-mild_day1"="lightblue", "AT2 cells-mild_day1"="lightblue", "AM (PBS)-mild_day1"="lightblue", 
                                                                                                 "B-lymphocytes-mild_day1"="lightblue", "T-lymphocytes-mild_day1"="lightblue", "MDM-mild_day1"="lightblue",
                                                                                                 "AT1 cells-mock"="lightgrey", "AT2 cells-mock"="lightgrey", "AM (PBS)-mock"="lightgrey", 
                                                                                                 "B-lymphocytes-mock"="lightgrey", "T-lymphocytes-mock"="lightgrey", "MDM-mock"="lightgrey"))
dev.off()

pdf("Ppara_target_genes_plot.pdf", width=13, height=13.5, onefile=F)
RidgePlot(merged_obj,features=c('Ppara3'), group.by='celltype_stage', fill.by="sample", cols = c("AT1 cells-severe_day1"="#CD5C5C", "AT2 cells-severe_day1"="#CD5C5C", "AM (PBS)-severe_day1"="#CD5C5C", 
                                                                                                 "B-lymphocytes-severe_day1"="#CD5C5C", "T-lymphocytes-severe_day1"="#CD5C5C", "MDM-severe_day1"="#CD5C5C",
                                                                                                 "AT1 cells-mild_day1"="lightblue", "AT2 cells-mild_day1"="lightblue", "AM (PBS)-mild_day1"="lightblue", 
                                                                                                 "B-lymphocytes-mild_day1"="lightblue", "T-lymphocytes-mild_day1"="lightblue", "MDM-mild_day1"="lightblue",
                                                                                                 "AT1 cells-mock"="lightgrey", "AT2 cells-mock"="lightgrey", "AM (PBS)-mock"="lightgrey", 
                                                                                                 "B-lymphocytes-mock"="lightgrey", "T-lymphocytes-mock"="lightgrey", "MDM-mock"="lightgrey"))
dev.off()

pdf("Pparg_target_genes_plot.pdf", width=13, height=13.5, onefile=F)
RidgePlot(merged_obj,features=c('Pparg4'), group.by='celltype_stage', fill.by="sample", cols = c("AT1 cells-severe_day1"="#CD5C5C", "AT2 cells-severe_day1"="#CD5C5C", "AM (PBS)-severe_day1"="#CD5C5C", 
                                                                                                 "B-lymphocytes-severe_day1"="#CD5C5C", "T-lymphocytes-severe_day1"="#CD5C5C", "MDM-severe_day1"="#CD5C5C",
                                                                                                 "AT1 cells-mild_day1"="lightblue", "AT2 cells-mild_day1"="lightblue", "AM (PBS)-mild_day1"="lightblue", 
                                                                                                 "B-lymphocytes-mild_day1"="lightblue", "T-lymphocytes-mild_day1"="lightblue", "MDM-mild_day1"="lightblue",
                                                                                                 "AT1 cells-mock"="lightgrey", "AT2 cells-mock"="lightgrey", "AM (PBS)-mock"="lightgrey", 
                                                                                                 "B-lymphocytes-mock"="lightgrey", "T-lymphocytes-mock"="lightgrey", "MDM-mock"="lightgrey"))
dev.off()


pdf("Ppard_target_genes_plot.pdf", width=13, height=13.5, onefile=F)
RidgePlot(merged_obj,features=c('Ppard5'), group.by='celltype_stage', fill.by="sample", cols = c("AT1 cells-severe_day1"="#CD5C5C", "AT2 cells-severe_day1"="#CD5C5C", "AM (PBS)-severe_day1"="#CD5C5C", 
                                                                                                 "B-lymphocytes-severe_day1"="#CD5C5C", "T-lymphocytes-severe_day1"="#CD5C5C", "MDM-severe_day1"="#CD5C5C",
                                                                                                 "AT1 cells-mild_day1"="lightblue", "AT2 cells-mild_day1"="lightblue", "AM (PBS)-mild_day1"="lightblue", 
                                                                                                 "B-lymphocytes-mild_day1"="lightblue", "T-lymphocytes-mild_day1"="lightblue", "MDM-mild_day1"="lightblue",
                                                                                                 "AT1 cells-mock"="lightgrey", "AT2 cells-mock"="lightgrey", "AM (PBS)-mock"="lightgrey", 
                                                                                                 "B-lymphocytes-mock"="lightgrey", "T-lymphocytes-mock"="lightgrey", "MDM-mock"="lightgrey"))
dev.off()
