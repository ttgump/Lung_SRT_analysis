rm(list=ls())
library(Seurat)
library(SeuratExtend)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Polychrome)
library(iTALK)
library(data.table)
library(RColorBrewer)

samples <- c("mock", "qing_day1", "qing_day2", "qing_day3", "qing_day5", "qing_day7", "qing_day14",
             "zhong_day1", "zhong_day2", "zhong_day3", "zhong_day5", "zhong_day7")
condition <- c("mock", "infected", "infected", "infected", "infected", "infected", "infected", "infected", "infected", "infected", "infected", "infected")

### load 16um data and merge 12 samples
object_list <- list(12)
for(i in 1:12) {
  object <- readRDS(paste("16um_SCT/",samples[i],"_SCT.rds",sep=""))
  DefaultAssay(object) <- "SCT"
  
  load(paste("16um_RCTD_deconvolution/", samples[i], "_myRCTD.RData", sep=""))
  
  results_df <- myRCTD@results$results_df
  results_df <- results_df[results_df$spot_class != "reject", ]
  results_df <- results_df[match(colnames(object), rownames(results_df)), ]
  rownames(results_df) <- colnames(object)
  results_df$first_type <- as.character(results_df$first_type)
  results_df$first_type[results_df$first_type %in% c("Cd163- Cd11c+ IMs", "Cd163+ Cd11c+ IMs", "AM (PBS)")] <- "AMs + IMs"
  results_df$first_type[results_df$first_type %in% c("Recruited macrophages", "M2 macrophages", "DCs", "Cd103+ DCs", "Ccl17+ DCs")] <- "MDM"
  
  object$condition <- samples[i]
  object$condition <- condition[i]
  
  cell_types <- factor(results_df$first_type)
  names(cell_types) <- colnames(object)
  object$cell_types <- cell_types
  
  object <- object[, object$cell_types %in% c("AT1 cells", "AT2 cells", "ABMC", "Activated AT2 cells", "AMs + IMs", "MDM", "T-lymphocytes", "Neutrophils")]
  
  object_list[[i]] <- object
}

names(object_list) <- samples
merged_object <- merge(object_list[[1]], object_list[2:12])
save(merged_object, file="iTALK_analysis_infected_vs_mock/16um_merged_selected_celltypes.RData")

# load("iTALK_analysis_infected_vs_mock/16um_merged_selected_celltypes.RData")

set.seed(42)
merged_object_sample <- merged_object

exp <- as.data.frame(t(as.matrix(merged_object_sample@assays$SCT$data)))
exp2 <- exp[, colSums(exp>0) > 0.005*nrow(exp)]
exp2$cell_type <- merged_object_sample$cell_types
exp2$compare_group <- merged_object_sample$condition
length(unique(merged_object_sample$cell_types))

### perform DE analysis in each cell types
deg_1 <-DEG(exp2 %>% filter(cell_type=='AT1 cells'),method='Wilcox',min_gene_expressed=0.01*nrow(exp2 %>% filter(cell_type=='AT1 cells')),min_valid_cells=10,contrast=c("infected","mock"))
deg_2 <-DEG(exp2 %>% filter(cell_type=='AT2 cells'),method='Wilcox',min_gene_expressed=0.01*nrow(exp2 %>% filter(cell_type=='AT2 cells')),min_valid_cells=10,contrast=c("infected","mock"))
deg_3 <-DEG(exp2 %>% filter(cell_type=='Activated AT2 cells'),method='Wilcox',min_gene_expressed=0.01*nrow(exp2 %>% filter(cell_type=='Activated AT2 cells')),min_valid_cells=10,contrast=c("infected","mock"))
deg_4 <-DEG(exp2 %>% filter(cell_type=='ABMC'),method='Wilcox',min_gene_expressed=0.01*nrow(exp2 %>% filter(cell_type=='ABMC')),min_valid_cells=10,contrast=c("infected","mock"))
deg_5 <-DEG(exp2 %>% filter(cell_type=='AMs + IMs'),method='Wilcox',min_gene_expressed=0.01*nrow(exp2 %>% filter(cell_type=='AMs + IMs')),min_valid_cells=10,contrast=c("infected","mock"))
deg_6 <-DEG(exp2 %>% filter(cell_type=='MDM'),method='Wilcox',min_gene_expressed=0.01*nrow(exp2 %>% filter(cell_type=='MDM')),min_valid_cells=10,contrast=c("infected","mock"))
deg_7 <-DEG(exp2 %>% filter(cell_type=='T-lymphocytes'),method='Wilcox',min_gene_expressed=0.01*nrow(exp2 %>% filter(cell_type=='T-lymphocytes')),min_valid_cells=10,contrast=c("infected","mock"))
deg_8 <-DEG(exp2 %>% filter(cell_type=='Neutrophils'),method='Wilcox',min_gene_expressed=0.01*nrow(exp2 %>% filter(cell_type=='Neutrophils')),min_valid_cells=10,contrast=c("infected","mock"))

### remove na and infinte values
deg_1 <- deg_1[!is.na(deg_1$q.value), ]
deg_2 <- deg_2[!is.na(deg_2$q.value), ]
deg_3 <- deg_3[!is.na(deg_3$q.value), ]
deg_4 <- deg_4[!is.na(deg_4$q.value), ]
deg_5 <- deg_5[!is.na(deg_5$q.value), ]
deg_6 <- deg_6[!is.na(deg_6$q.value), ]
deg_7 <- deg_7[!is.na(deg_7$q.value), ]
deg_8 <- deg_8[!is.na(deg_8$q.value), ]

deg_1$logFC[deg_1$logFC==Inf] <- 10
deg_1$logFC[deg_1$logFC==-Inf] <- -10
deg_2$logFC[deg_2$logFC==Inf] <- 10
deg_2$logFC[deg_2$logFC==-Inf] <- -10
deg_3$logFC[deg_3$logFC==Inf] <- 10
deg_3$logFC[deg_3$logFC==-Inf] <- -10
deg_4$logFC[deg_4$logFC==Inf] <- 10
deg_4$logFC[deg_4$logFC==-Inf] <- -10
deg_5$logFC[deg_5$logFC==Inf] <- 10
deg_5$logFC[deg_5$logFC==-Inf] <- -10
deg_6$logFC[deg_6$logFC==Inf] <- 10
deg_6$logFC[deg_6$logFC==-Inf] <- -10
deg_7$logFC[deg_7$logFC==Inf] <- 10
deg_7$logFC[deg_7$logFC==-Inf] <- -10
deg_8$logFC[deg_8$logFC==Inf] <- 10
deg_8$logFC[deg_8$logFC==-Inf] <- -10

save(list=c("deg_1", "deg_2", "deg_3", "deg_4", "deg_5", "deg_6", "deg_7", "deg_8"), file="iTALK_analysis_infected_vs_mock/16um_iTalk_wilcox_DE_8celltypes.RData")


rownames(deg_1) <- toupper(rownames(deg_1))
deg_1$gene <- toupper(deg_1$gene)
rownames(deg_2) <- toupper(rownames(deg_2))
deg_2$gene <- toupper(deg_2$gene)
rownames(deg_3) <- toupper(rownames(deg_3))
deg_3$gene <- toupper(deg_3$gene)
rownames(deg_4) <- toupper(rownames(deg_4))
deg_4$gene <- toupper(deg_4$gene)
rownames(deg_5) <- toupper(rownames(deg_5))
deg_5$gene <- toupper(deg_5$gene)
rownames(deg_6) <- toupper(rownames(deg_6))
deg_6$gene <- toupper(deg_6$gene)
rownames(deg_7) <- toupper(rownames(deg_7))
deg_7$gene <- toupper(deg_7$gene)
rownames(deg_8) <- toupper(rownames(deg_8))
deg_8$gene <- toupper(deg_8$gene)

comm_list<-c('growth factor','other','cytokine','checkpoint')

cell_cols<-as.character(brewer.pal(8, "Set3"))
names(cell_cols) <- c("AT1 cells", "AT2 cells", "ABMC", "Activated AT2 cells", "AMs + IMs", "MDM", "T-lymphocytes", "Neutrophils")


### combine all DE analysis and plot
all.res <- NULL
for(comm_type in comm_list){
  res<-NULL
  sig.res<-NULL
  pdf(paste(comm_type, "_16um_iTalk.pdf", sep=""), width=15, height=8, onefile=F)
  par(mfrow=c(1,2))
  
  for(i in 1:7){
    for(j in (i+1):8) {
      res_cat<-FindLR(get(paste("deg_", i, sep="")),get(paste("deg_", j, sep="")),datatype='DEG',comm_type=comm_type)
      res_cat<-res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC,decreasing=T),]
      res_cat<-res_cat[sign(res_cat$cell_from_logFC)==sign(res_cat$cell_to_logFC), ]
      if(nrow(res_cat) > 0) {
        res_cat$GainOrLoss <- "red"
        res_cat$GainOrLoss[res_cat$cell_from_logFC < 0] <- "blue"
        res<-rbind(res,res_cat)
      }
    }
  }
  
  if(!is.null(res))
  {
    sig.res<-res[order(res$cell_from_logFC,decreasing=T),]
    sig.res<-sig.res[1:min(50, nrow(sig.res)),]
    
    NetView(res,col=cell_cols,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    LRPlot(sig.res,datatype='DEG',
           cell_col=cell_cols,link.arr.lwd=abs(sig.res$cell_from_logFC),
           link.arr.width=abs(sig.res$cell_from_logFC), 
           link.arr.col=sig.res$GainOrLoss)
    title(comm_type)
  }
  
  
  dev.off()
  #  write.csv(res, file=paste(comm_type, "_covid19_iTalk.csv", sep=""), row.names=F)
  all.res <- rbind(all.res, res)
}


sig.res<-all.res[order(all.res$cell_from_logFC,decreasing=T),]
sig.res<-sig.res[1:min(80, nrow(sig.res)),]

pdf("all_pathway_16um_iTalk.pdf", width=15, height=8, onefile=F)
par(mfrow=c(1,2))
NetView(all.res,col=cell_cols,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
LRPlot(sig.res,datatype='DEG',
       cell_col=cell_cols,link.arr.lwd=abs(sig.res$cell_from_logFC),
       link.arr.width=abs(sig.res$cell_from_logFC), 
       link.arr.col=sig.res$GainOrLoss)
title("all pathways")
dev.off()

write.csv(all.res, file="all_pathway_16um_iTalk.csv", row.names=F)
