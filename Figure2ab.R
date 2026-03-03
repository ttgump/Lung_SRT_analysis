rm(list=ls())
library(Seurat)
library(semla)
library(tibble)
library(ggplot2)
library(patchwork)
library(scico)
library(tidyr)
library(dplyr)
library(viridis)

samples <- c("mock", "mild_day1", "mild_day2", "mild_day3", "mild_day5", "mild_day7", "mild_day14",
             "severe_day1", "severe_day2", "severe_day3", "severe_day5", "severe_day7")
dirs <- c("mock/outs/binned_outputs/", 
          "mild_day1/outs/binned_outputs/",
          "mild_day2/outs/binned_outputs/", 
          "mild_day3/outs/binned_outputs/",
          "mild_day5/outs/binned_outputs/", 
          "mild_day7/outs/binned_outputs/",
          "mild_day14/outs/binned_outputs/", 
          "severe_day1/outs/binned_outputs/",
          "severe_day2/outs/binned_outputs/", 
          "severe_day3/outs/binned_outputs/",
          "severe_day5/outs/binned_outputs/", 
          "severe_day7/outs/binned_outputs/")

##### Custom Functions ####
RunGeneDistanceCor <- function(
    object, 
    dist_colname,
    genes_select, 
    dist_cutoff = NULL,
    dist_lower_cutoff = 0){
  
  if (is.null(dist_cutoff)) dist_cutoff <- 500
  
  gene_data <- object[[]] %>% 
    bind_cols(FetchData(object, vars = genes_select)) %>% 
    filter(.[[dist_colname]] < dist_cutoff & .[[dist_colname]] >= dist_lower_cutoff)
  
  message(paste0("Selecting ", length(genes_select), " genes"))
  message(paste0("Selecting ", nrow(gene_data), " spots"))
  
  message()
  
  gene_cor <- setNames(lapply(genes_select, function(g){
    cor_res <- cor.test(x = gene_data[[dist_colname]], y = gene_data[[g]])
    cor_res <- cbind(cor_res$estimate, cor_res$p.value)
    colnames(cor_res) <- c("cor", "pval")
    rownames(cor_res) <- g
    return(cor_res)
  }), nm = genes_select)
  
  gene_cor_df <- do.call("rbind", gene_cor) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene") %>% 
    as_tibble() %>%
    filter(!is.na(cor)) %>% 
    arrange(-cor)
  gene_cor_df$pval_FDR <- p.adjust(gene_cor_df$pval, method = "BH")
  gene_cor_df$sign_005 <- ifelse(gene_cor_df$pval < 0.05, TRUE, FALSE)
  gene_cor_df$sign_001 <- ifelse(gene_cor_df$pval < 0.01, TRUE, FALSE)
  gene_cor_df$sign_padj_005 <- ifelse(gene_cor_df$pval_FDR < 0.05, TRUE, FALSE)
  gene_cor_df$sign_padj_001 <- ifelse(gene_cor_df$pval_FDR < 0.01, TRUE, FALSE)
  
  return(gene_cor_df)
}


#### custom colors
col_scale_mako <- viridis::mako(13, direction = -1)[1:10]
col_scale_acton <- scico::scico(10, direction = -1, palette = "acton")
col_scale_div_custom <- c(rev(col_scale_mako[1:8]), "white", col_scale_acton[1:8])

load("../mock_mild_7_samples_ATs.RData")
load("../mock_severe_6_samples_ATs.RData")

se.hd.list <- list(0)

j <- 1
for(i in c(1:12)) {
  st.dir <- file.path(dirs[i])
  res <- "08um"
  
  res.dir <- list.dirs(st.dir, recursive = FALSE) |> 
    stringr::str_subset(pattern = res)
  
  infoTable <- data.frame(
    samples = list.files(res.dir,
                         full.names = TRUE, recursive = TRUE,
                         pattern = paste0("^filtered_feature.+.h5$")
    ),
    spotfiles = list.files(res.dir,
                           full.names = TRUE, recursive = TRUE,
                           pattern = "parquet$|positions.csv$"
    ),
    imgs = list.files(res.dir,
                      recursive = TRUE,
                      full.names = TRUE, pattern = "hires"
    ),
    json = list.files(res.dir,
                      recursive = TRUE,
                      full.names = TRUE, pattern = "^scalefactors"
    ),
    resolution = res,
    sample_ID = samples[i]
  )
  
  se.hd <- ReadVisiumData(infoTable)
  se.hd <- LoadImages(se.hd)
  
  load(paste("../", samples[i], "_RCTD_Plots/results.RData", sep=""))
  
  results_df <- results$results_df
  results_df <- results_df[results_df$spot_class != "reject", ]
  results_df <- results_df[match(colnames(se.hd), rownames(results_df)), ]
  rownames(results_df) <- colnames(se.hd)
  
  if(i <= 7) {
    current_obj <- mild_obj[, mild_obj$sample == samples[i]]
  } else {
    current_obj <- severe_obj[, severe_obj$sample == samples[i]]
  }
  
  current_obj_cell_type <- current_obj$cell_type2
  names(current_obj_cell_type) <- substr(colnames(current_obj), nchar(samples[i])+2, nchar(colnames(current_obj)))
  
  results_df$first_type <- as.character(results_df$first_type)
  results_df[names(current_obj_cell_type), "first_type"] <- current_obj_cell_type
  
  cell_type <- as.character(results_df$first_type)
  cell_type[cell_type %in% c("Cd163- Cd11c+ IMs", "Cd163+ Cd11c+ IMs", "AM (PBS)")] <- "AMs + IMs"
  cell_type[cell_type %in% c("Recruited macrophages", "M2 macrophages", "DCs", "Cd103+ DCs", "Ccl17+ DCs")] <- "MDM"
  names(cell_type) <- rownames(results_df)
  cell_type <- factor(cell_type)
  
  ### if there is no ABMC in current sample, jump to the next sample
  if(sum(cell_type=="ABMC", na.rm = T) == 0) {
    next
  }
  
  se.hd$cell_type <- cell_type
  
  se.hd <- RadialDistance(se.hd, column_name = "cell_type", selected_groups = "ABMC", convert_to_microns = TRUE)
  if(sum(is.na(se.hd$`r_dist_ABMC`)) >= 100) {
    next
  }
  
  j <- j+1
  se.hd.list[[j]] <- se.hd
  
  ### plot cell type distributions around ABMC cells
  rdist_data <- se.hd@meta.data
  rdist_data <- rdist_data |> filter(!is.na(`r_dist_ABMC`))
  
  rdist_data$`r_dist_ABMC` |> max()
  rdist_data$`r_dist_ABMC` |> min()
  rdist_data <- rdist_data |> filter(cell_type %in% c("AT1 cells", "AT2 cells", "Activated AT2 cells", "ABMC"))
  
  min.dist <- min(rdist_data$`r_dist_ABMC`)
  max.dist <- max(rdist_data$`r_dist_ABMC`)
  bin_breaks <- seq(min.dist, max.dist, 100)
  
  rdist_data$dist_bin <- cut(rdist_data$`r_dist_ABMC`, bin_breaks, 
                                  labels = bin_breaks[-length(bin_breaks)]+50)
  
  rdist_data$dist_bin <- rdist_data$dist_bin |> as.character() |> as.numeric()
  
  rdist_data$dist_bin_min <- rdist_data$dist_bin - 50
  rdist_data$dist_bin_max <- rdist_data$dist_bin + 50
  
  dist_cutoff_plot <- 2000
  
  rdist_data_filter <- rdist_data |> filter(`r_dist_ABMC` < dist_cutoff_plot)
  rdist_data_filter_1 <- rdist_data_filter |> filter(cell_type %in% c("AT1 cells", "AT2 cells", "Activated AT2 cells", "ABMC"))
  rdist_data_filter_1$cell_type <- factor(rdist_data_filter_1$cell_type, levels = c("ABMC", "Activated AT2 cells", "AT2 cells", "AT1 cells"))
  
  pdf(paste(samples[i], "_celltypes.pdf", sep=""), width=4.5, height=5, onefile=F)  # !Main Figure
  print(ggplot(rdist_data_filter_1, aes(x=`r_dist_ABMC`)) +
    geom_vline(xintercept = 0, color="orange") +
    geom_density(n=40, bounds=c(0, Inf)) +
    facet_wrap(~cell_type, nrow = 6, scales = "free_y") +
    theme_classic() +
    theme(strip.text.y.right = element_text(angle=0, hjust=0), 
          strip.background = element_rect(color = NA), 
          legend.position = "none",
          axis.text.x = element_text(angle=45, vjust = 1, hjust=1), 
          # axis.title.y = element_blank(),
          panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
    ))
  dev.off()
  

  ### plot gene distributions around ABMC cells
  se.hd <-  se.hd |>
    NormalizeData() |>
    ScaleData() |>
    FindVariableFeatures()
  
  gene_selected <- c(VariableFeatures(se.hd), "Cxcl10", "Ccl20", "Cxcl5", "Ifit1bl1", "Cxcl10", "Hp", "Rsad2", "Cxcl9", "Ifi205",
                     "Mxd1", "Icam1", "Cxcl5", "Isg15", "Ccl2", "Ccl3")
  gene_selected <- unique(gene_selected)
  gene_cor_df <- RunGeneDistanceCor(object = se.hd, dist_colname = "r_dist_ABMC", genes_select = gene_selected, dist_cutoff = 2000)
  
  dist_cutoff_plot <- 2000
  
  # Plot top n positive correlated genes
  n_genes <- 20
  genes_pos <- gene_cor_df %>% arrange(desc(cor)) %>% head(n_genes) %>% filter(sign_001 == TRUE) %>% pull(gene)
  
  p_dat_g_pos_line <- se.hd[[]] %>% 
    bind_cols(FetchData(se.hd, vars = c(genes_pos))) %>% 
    filter(`r_dist_ABMC` < dist_cutoff_plot) %>% 
    pivot_longer(all_of(c(genes_pos)), names_to = "gene", values_to = "expression") %>% 
    mutate_at(.vars="gene", ~ factor(., levels = c(genes_pos)))
  p_dat_g_pos_line$r_dist <- round(p_dat_g_pos_line$`r_dist_ABMC`)
  
  p_dat_g_pos_cor <- gene_cor_df %>% arrange(desc(cor)) %>% filter(gene %in% genes_pos)
  p_dat_g_pos_cor$gene <- factor(p_dat_g_pos_cor$gene, levels = genes_pos)
  
  # Plot top n negative correlated genes
  n_genes <- 20
  genes_neg_df1 <- gene_cor_df %>% arrange(cor) %>% head(n_genes) %>% filter(sign_001 == TRUE)
  
  genes_neg_df2 <- gene_cor_df[gene_cor_df$gene %in% c("Cxcl10", "Ccl20", "Cxcl5", "Ifit1bl1", "Cxcl10", "Hp", "Rsad2", "Cxcl9", "Ifi205",
                     "Mxd1", "Icam1", "Cxcl5", "Isg15", "Ccl2", "Ccl3"), ] %>% arrange(cor) %>% filter(sign_001 == TRUE)
  genes_neg_df <- rbind(genes_neg_df1, genes_neg_df2)
  genes_neg_df <- genes_neg_df[genes_neg_df$cor < 0, ]
  genes_neg_df <- genes_neg_df[!duplicated(genes_neg_df), ]
  genes_neg_df <- genes_neg_df[order(genes_neg_df$cor),]
  genes_neg <- genes_neg_df %>% pull(gene)
  
  p_dat_g_neg_line <- se.hd[[]] %>% 
    bind_cols(FetchData(se.hd, vars = c(genes_neg))) %>% 
    filter(`r_dist_ABMC` < dist_cutoff_plot) %>% 
    pivot_longer(all_of(c(genes_neg)), names_to = "gene", values_to = "expression") %>% 
    mutate_at(.vars="gene", ~ factor(., levels = c(genes_neg)))
  p_dat_g_neg_line$r_dist <- round(p_dat_g_neg_line$`r_dist_ABMC`)
  
  p_dat_g_neg_cor <- gene_cor_df %>% arrange((cor)) %>% filter(gene %in% genes_neg)
  p_dat_g_neg_cor$gene <- factor(p_dat_g_neg_cor$gene, levels = genes_neg)
  
  # plot
  # pos:
  p_g_pos_line <- ggplot(p_dat_g_pos_line, aes(`r_dist_ABMC`, expression)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth=0.5, fill="grey80", color = "black") +
    geom_vline(aes(xintercept = 0, color = "border"), 
               linetype = "dashed", linewidth=0.5, color="black") +
    facet_wrap(~gene, scales = "free_y", ncol = 1, strip.position = "right") +
    theme_classic() +
    theme(strip.text.y.right = element_text(angle=0, hjust=0), 
          strip.background = element_rect(color = NA), 
          legend.position = "none",
          axis.text.x = element_text(angle=45, vjust = 1, hjust=1), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          # axis.title.y = element_blank(),
          panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
    ); p_g_pos_line
  
  cor_lims <- 0.51
  p_g_pos_cor <- ggplot(p_dat_g_pos_cor, aes(x="cor", y=reorder(gene, cor), fill=cor)) +
    geom_tile(height = 0.75) +
    scale_fill_gradientn(colours = col_scale_div_custom, limits = c(-cor_lims, cor_lims)) +
    theme_classic() +
    theme(text = element_text(color="black", size=10),
          axis.title = element_blank(), 
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "left",
          legend.text = element_text(color="black", size=10)); p_g_pos_cor
  
  p1 <- (p_g_pos_cor | p_g_pos_line) + plot_layout(widths = c(1,3))
  
  # neg:
  p_g_neg_line <- ggplot(p_dat_g_neg_line, aes(`r_dist_ABMC`, expression)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth=0.5, fill="grey80", color = "black") +
    geom_vline(aes(xintercept = 0, color = "border"), 
               linetype = "dashed", linewidth=0.5, color="black") +
    facet_wrap(~gene, scales = "free_y", ncol = 1, strip.position = "right") +
    theme_classic() +
    theme(strip.text.y.right = element_text(angle=0, hjust=0), 
          strip.background = element_rect(color = NA), 
          legend.position = "none",
          axis.text.x = element_text(angle=45, vjust = 1, hjust=1), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          # axis.title.y = element_blank(),
          panel.grid.major.x = element_line(linewidth=0.5, color="grey90")
    ); p_g_neg_line
  
  cor_lims <- 0.51
  p_g_neg_cor <- ggplot(p_dat_g_neg_cor, aes(x="cor", y=reorder(gene, desc(cor)), fill=cor)) +
    geom_tile(height = 0.75) +
    scale_fill_gradientn(colours = col_scale_div_custom, limits = c(-cor_lims, cor_lims)) +
    theme_classic() +
    theme(text = element_text(color="black", size=10),
          axis.title = element_blank(), 
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "left",
          legend.text = element_text(color="black", size=10)); p_g_neg_cor
  
  p2 <- (p_g_neg_cor | p_g_neg_line) + plot_layout(widths = c(1,3))
  
  
  # Export
  pdf(paste(samples[i], "_ABMC_dist_gene_pos.pdf", sep=""), width=4, height=3+(length(genes_pos)-5)*2.5/15, onefile=F)  # !Main Figure
  print(p1)
  dev.off()
  
  pdf(paste(samples[i], "_ABMC_dist_gene_neg.pdf", sep=""), width=4, height=3+(length(genes_neg)-5)*2.5/15, onefile=F)  # !Main Figure
  print(p2)
  dev.off()
}

