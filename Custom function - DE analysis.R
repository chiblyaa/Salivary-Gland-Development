# This script will perform analysis associated with Figures 5-7
# First we'll perform integratio of all datasets for visualization
# and trajectory analysis. Then, we perform differential expression
# analysis between specific populations and cross-reference with a 
# database of transcription factors. Lastly, we perform correlation
# analysis for Bhlha15 to identify genes associated with acinar lineage.

## Custom function to perform differential expression analysis between populations
transition.analysis <- function(seuratobj, ident1, ident2, tfs, cols, filename, up=5, down=5, group.by="stage"){
  ## The use of this function requires the following criteria to work properly:
  ## seuratobj: annotated seurat object labeled with cell identities to compare and with metadata for 'stage'
  ## ident1 and ident2: identities to compare.
  ## tfs: a table of transcription factors with a 'Gene' column. Alternatively use list of genes of interest
  ## cols: a combination of TWO colors for violin plots
  ## filename: a desired name for output files. The function will create list of genes and violin plot
  ## up and down: number of genes to show in violin plot
  require(Seurat)
  require(cowplot)
  require(dplyr)
  require(ggplot2)
  DefaultAssay(seuratobj) <- "RNA" #For differential expression we must use the raw non-integrated data
  comparison <- FindMarkers(object = seuratobj, ident.1 = ident1, ident.2 = ident2, only.pos = F)
  comparison <- comparison[comparison$p_val_adj<0.05, ]
  write.csv(comparison, file = paste0(filename, ".csv"), row.names = T)
  comparison.tfs <- comparison[rownames(comparison) %in% tfs$Gene, ]
  write.csv(comparison.tfs, file = paste0(filename, "-tfs.csv"), row.names = T)
  comparison$gene <- rownames(comparison)
  comparison.tfs$gene <- rownames(comparison.tfs)
  toptfs <- rbind((comparison.tfs %>% top_n(up,avg_logFC)), (comparison.tfs %>% top_n(-down,avg_logFC)))
  toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
  plot.list <- list()
  if(length(toptfs$gene!=0)){
  for(i in 1:length(toptfs$gene)){
    plot.list[[i]] <- VlnPlot(seuratobj, idents = c(ident1, ident2), cols = cols, features = toptfs$gene[i], pt.size = 0, group.by = group.by) + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 12, face = "italic", hjust = 0.07, vjust = -2), plot.margin = margin(unit = "pt", t=0)) +NoLegend()
  }
  pdf(paste0(filename, ".pdf"), width = 5, height = 2,useDingbats = F) 
  print(plot_grid(plotlist = plot.list, ncol = 5))
  dev.off() 
  }
  results <- list(comparison, comparison.tfs, toptfs)
  return(results)
}
