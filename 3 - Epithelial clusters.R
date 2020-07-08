library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(superheat)
library(Hmisc)


########################## Load previously annotated SEURAT files############################
#############################################################################################

e12smg <- readRDS(file = "~/Import to SEURAT/E12 SMG Annotated (SEURAT v3).rds")
e14smg <- readRDS(file = "~/Import to SEURAT/E14 SMG Annotated (SEURAT v3).rds")
e16smg <- readRDS(file = "~/Import to SEURAT/E16 SMG Annotated (SEURAT v3).rds")
p1smg <- readRDS(file = "~/Import to SEURAT/p1 SMG Annotated (SEURAT v3).rds")
adultsmg <- readRDS(file = "~/Adult SMG (ICR CH3)/Adult SMG (ICR-CH3) Annotated (SEURAT v3).rds")

##### Verify consistency in Cell Type annotation - correct if necessary

DimPlot(e12smg, pt.size = 1, group.by = "CellType", label = T,label.size = 6)
DimPlot(e14smg, pt.size = 1, group.by = "CellType", label = T,label.size = 6)
DimPlot(e16smg, pt.size = 1, group.by = "CellType", label = T,label.size = 6)
DimPlot(p1smg, pt.size = 1, group.by = "CellType", label = T,label.size = 6)
DimPlot(adultsmg, pt.size = 1, group.by = "CellType", label = T,label.size = 6)

##### Add stage metadata if it hasn't been added before

e12smg[["stage"]] <- "E12"
e14smg[["stage"]] <- "E14"
e16smg[["stage"]] <- "E16"
p1smg[["stage"]] <- "P1"
adultsmg[["stage"]] <- "Adult"

#### Separate epithelial clusters using Cell Type ID

e12smg <- SetIdent(e12smg, value = "CellType")
e14smg <- SetIdent(e14smg, value = "CellType")
e16smg <- SetIdent(e16smg, value = "CellType")
p1smg <- SetIdent(p1smg, value = "CellType")
adultsmg <- SetIdent(adultsmg, value = "CellType")

e12epithelium <- subset(e12smg,idents = c("End bud", "Krt19+ duct"))
DimPlot(e12epithelium)

e14epithelium <- subset(e14smg,idents = c("End bud", "Krt19+ duct"))
DimPlot(e14epithelium)

e16epithelium <- subset(e16smg,idents = c("End bud", "Krt19+ duct", "Basal duct", "Myoepithelial"))
DimPlot(e16epithelium)

### Rename Mitotic cells to SMGC+ in P1 dataset for consistency
p1smg <- RenameIdents(p1smg, "Mitotic cells" = "Smgc+ Proacinar")
p1epithelium <- subset(p1smg,idents = c("Bpifa2+ Proacinar","Smgc+ Proacinar", "Krt19+ duct", "Basal duct", "Myoepithelial"))
DimPlot(p1epithelium)

adultepithelium <- subset(adultsmg,idents = c("Smgc+","Acinar","Bpifa2+ Acinar", "Striated duct", "Basal duct", "Myoepithelial", "Ascl3", "Intercalated duct", "GCT"))
DimPlot(adultepithelium)

#### Clear unnecessary files to save memory
remove(e12smg, e14smg, e16smg, p1smg, adultsmg)

#############################################################################################
########################## Re-cluster epithelium with SEURAT ################################

###############################################################
########################## E12 ################################

e12epithelium <- NormalizeData(object = e12epithelium, normalization.method = "LogNormalize", scale.factor = 10000)
e12epithelium <- FindVariableFeatures(object = e12epithelium, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = e12epithelium), 10)
plot1 <- VariableFeaturePlot(object = e12epithelium)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(x = e12epithelium)
e12epithelium <- ScaleData(object = e12epithelium, features = all.genes, vars.to.regress = "percent.mt")
e12epithelium <- RunPCA(object = e12epithelium, features = VariableFeatures(object = e12epithelium))
ElbowPlot(object = e12epithelium)
e12epithelium <- FindNeighbors(object = e12epithelium, dims = 1:10,force.recalc = T)
e12epithelium <- FindClusters(object = e12epithelium, resolution = 0.5)
e12epithelium <- RunTSNE(object = e12epithelium, dims = 1:10)

DimPlot(object = e12epithelium, reduction = "tsne", pt.size = 1.5,label = T,label.size = 8, cols = c("skyblue1", "red", "darkred", "skyblue3", "royalblue2", "mediumpurple1", "hotpink", "gray")) + theme(axis.text = element_text(size = 14, colour = "black"))

### Identify E12 epithelium Unsupervised markers
e12epithelium <- SetIdent(e12epithelium, value = "seurat_clusters")
e12epi.markers <- FindAllMarkers(e12epithelium, logfc.threshold = 0.25, only.pos = T)
e12epi.markers.top5 <- e12epi.markers %>% group_by(cluster) %>% top_n(-5,p_val_adj)
DotPlot(object = e12epithelium, features = c(rev(unique(e12epi.markers.top5$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(angle = 90, size = 14))
write.table(e12epi.markers, file = "E12 Unsupervised Epithelium markers.txt", sep = "\t", row.names = F)

### Generate tSNE plots for Figure

plot1 <- FeaturePlot(e12epithelium, features = c("Krt5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(e12epithelium, features = c("Krt19"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(e12epithelium, features = c("Krt14"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot4 <- FeaturePlot(e12epithelium, features = c("Sox10"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes()+ theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3, plot4))

plot1 <- FeaturePlot(e12epithelium, features = c("Ung"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(e12epithelium, features = c("Aurka"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(e12epithelium, features = c("Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)

plot1 <- FeaturePlot(e12epithelium, features = c("Snhg4"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt")) 
plot2 <- FeaturePlot(e12epithelium, features = c("Krt15"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot3 <- FeaturePlot(e12epithelium, features = c("Wnt6"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot4 <- FeaturePlot(e12epithelium, features = c("Ccne1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot5 <- FeaturePlot(e12epithelium, features = c("Ccnb1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot6 <- FeaturePlot(e12epithelium, features = c("Col1a2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot7 <- FeaturePlot(e12epithelium, features = c("Eif4a1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7), ncol = 1)

FeaturePlot(e12epithelium, features = c("Trp63"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes()
FeaturePlot(e12epithelium, features = c("Aurka"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes()
FeaturePlot(e12epithelium, features = c("Ung"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes()
FeaturePlot(e12epithelium, features = c("Mki67"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes()
FeaturePlot(e12epithelium, features = c("Kit"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes()



### Annotate epithelial clusters
e12epithelium <- RenameIdents(e12epithelium, '1' = "Krt19+ duct",'0' = "End bud",'3'= "End bud",'2'= "Krt19+ duct",'4'= "End bud",'5'= "Mesenchyme",'6'= "Mitotic cells")
Idents(object = e12epithelium) <- factor(x = Idents(object = e12epithelium), levels = c("End bud", "Krt19+ duct", "Mesenchyme", "Mitotic cells"))
e12epithelium[["Epi.CellType"]] <- Idents(object = e12epithelium)
e12epithelium <- SetIdent(e12epithelium, value = "Epi.CellType")
DimPlot(object = e12epithelium, reduction = "tsne", pt.size = 1.5,label = T,label.size = 8, cols = c("skyblue3", "red", "mediumpurple1", "hotpink")) + theme(axis.text = element_text(size = 14, colour = "black"))

### Determine differentially expressed genes in duct vs end bud clusters
e12epithelium <- SetIdent(e12epithelium, value = "Epi.CellType")
e12.celltype.markers <- FindAllMarkers(e12epithelium,only.pos = T)
e12.celltype.markers <- e12.celltype.markers[e12.celltype.markers$p_val_adj<0.05,]
e12.celltype.markers$gene <- rownames(e12.celltype.markers)
write.table(e12.celltype.markers, file = "E12 Cell Type markers (SEURAT v3).txt", sep = "\t", row.names = F)
e12.celltype.markers.top10 <- e12.celltype.markers %>% group_by(cluster) %>% top_n(-10, p_val_adj)
DotPlot(object = e12epithelium, features = c(rev(e12.celltype.markers.top10$gene), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(size = 14))
## Belinda's markers
DotPlot(e12epithelium, group.by = "Epi.CellType", features = rev(c("BC006965", "Ngfr", "Tfap2b", "Shisa2", "Col9a1", "Etv4", "Hpgd", "Etv5", "Ccnd1", "Selm", "Krt15", "Krt5", "Upk3bl", "Igfbp5", "Krt19", "Anxa2", "Mafb", "Id1", "Wnt4", "Sfn")), cols = c("Spectral")) + theme(axis.text.x =  element_text(angle = 90))
saveRDS(e12epithelium, file = "E12 Epithelium annotated (SEURAT v3).rds")

###############################################################
########################## E14 ################################
#### Re-cluster epithelium with SEURAT 

e14epithelium <- NormalizeData(object = e14epithelium, normalization.method = "LogNormalize", scale.factor = 10000)
e14epithelium <- FindVariableFeatures(object = e14epithelium, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = e14epithelium), 10)
plot1 <- VariableFeaturePlot(object = e14epithelium)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(x = e14epithelium)
e14epithelium <- ScaleData(object = e14epithelium, features = all.genes, vars.to.regress = "percent.mt")
e14epithelium <- RunPCA(object = e14epithelium, features = VariableFeatures(object = e14epithelium))
ElbowPlot(object = e14epithelium)
e14epithelium <- FindNeighbors(object = e14epithelium, dims = 1:10, force.recalc = T)
e14epithelium <- FindClusters(object = e14epithelium, resolution = 0.6)
e14epithelium <- RunTSNE(object = e14epithelium, dims = 1:10)
DimPlot(object = e14epithelium, group.by = "seurat_clusters", reduction = "tsne", pt.size = 1,label = F,label.size = 6, cols = c("blue", "skyblue3", "royalblue2", "turquoise", "mediumpurple2", "orange", "darkred")) + NoLegend() +theme()
### Identify E14 epithelium Unsupervised markers
e14epithelium <- SetIdent(e14epithelium, value = "seurat_clusters")
e14epi.markers <- FindAllMarkers(object = e14epithelium, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.25)
e14epi.markers <- e14epi.markers[e14epi.markers$p_val_adj<0.05,]
write.table(e14epi.markers, file = "E14 Unsupervised Epithelium markers.txt", sep = "\t", row.names = F)
e14epi.markers.top5 <- e14epi.markers %>% group_by(cluster) %>% top_n(-5,p_val_adj)
DotPlot(object = e14epithelium, features = c(rev(unique(e14epi.markers.top5$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(angle = 90, size = 14))
### Generate tSNE plots for Figure
plot1 <- FeaturePlot(e14epithelium, features = c("Krt5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(e14epithelium, features = c("Krt7"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot3 <- FeaturePlot(e14epithelium, features = c("Krt14"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot4 <- FeaturePlot(e14epithelium, features = c("Sox10"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3, plot4))

plot1 <- FeaturePlot(e14epithelium, features = c("Ung"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(e14epithelium, features = c("Aurka"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(e14epithelium, features = c("Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)

plot1 <- FeaturePlot(e14epithelium, features = c("Actg1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(e14epithelium, features = c("Ccne1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot3 <- FeaturePlot(e14epithelium, features = c("Cenpf"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot4 <- FeaturePlot(e14epithelium, features = c("Ccnd1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot5 <- FeaturePlot(e14epithelium, features = c("Col1a2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot6 <- FeaturePlot(e14epithelium, features = c("Foxq1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot7 <- FeaturePlot(e14epithelium, features = c("Id3"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7), ncol = 1)

FeaturePlot(e14epithelium, features = c("Krt5"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))

## Annotate epithelial clusters
e14epithelium <- RenameIdents(e14epithelium, '0' = "End bud", '1' = "End bud",'2'= "End bud",'3'= "End bud",'4'= "Mesenchyme",'5'= "Duct",'6'= "Duct")
Idents(object = e14epithelium) <- factor(x = Idents(object = e14epithelium), levels = c("End bud", "Duct", "Mesenchyme"))
e14epithelium[["Epi.CellType"]] <- Idents(object = e14epithelium)
e14epithelium <- SetIdent(e14epithelium, value = "Epi.CellType")
DimPlot(object = e14epithelium, reduction = "tsne", pt.size = 1.5,label = T,label.size = 8, cols = c("skyblue3", "red", "mediumpurple1", "hotpink")) + theme(axis.text = element_text(size = 14, colour = "black"))
### Determine differentially expressed genes in duct vs end bud clusters

e14.celltype.markers <- FindAllMarkers(e14epithelium, only.pos = T)
e14.celltype.markers <- e14.celltype.markers[e14.celltype.markers$p_val_adj<0.05,]
e14.celltype.markers$gene <- rownames(e14.celltype.markers)
write.table(e14.celltype.markers, file = "E14 cell type markers (SEURAT v3).txt", sep = "\t", row.names = F)
e14.celltype.markers.top10 <- e14.celltype.markers %>% group_by(cluster) %>% top_n(-10, p_val_adj)
DotPlot(object = e14epithelium, features = c(rev(e14.celltype.markers.top10$gene), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(size = 14))
## Belinda's markers
DotPlot(e14epithelium, group.by = "Epi.CellType", features = rev(c("Sema3b", "Ccl20", "Mgp", "Col9a1", "Id4", "Dcpp1", "Tinagl1", "Myl9", "Etv5", "Arhgdib", "Foxq1", "Ppp1r1b", "Nupr1", "Krt7", "Cyp2f2", "Tacstd2", "Wfdc18", "Runx1", "Krt19", "Wnt5b")), cols = c("Spectral")) + theme(axis.text.x =  element_text(angle = 90))
saveRDS(e14epithelium, file = "E14 Epithelium annotated (SEURAT v3).rds")

######################################################################
########################## E16 #######################################
#############  Re-cluster epithelium with SEURAT 

e16epithelium <- NormalizeData(object = e16epithelium, normalization.method = "LogNormalize", scale.factor = 10000)
e16epithelium <- FindVariableFeatures(object = e16epithelium, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = e16epithelium), 10)
plot1 <- VariableFeaturePlot(object = e16epithelium)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(x = e16epithelium)
e16epithelium <- ScaleData(object = e16epithelium, features = all.genes)
e16epithelium <- RunPCA(object = e16epithelium, features = VariableFeatures(object = e16epithelium))
ElbowPlot(object = e16epithelium)
e16epithelium <- FindNeighbors(object = e16epithelium, dims = 1:10, force.recalc = T)
e16epithelium <- FindClusters(object = e16epithelium, resolution = 0.8)
e16epithelium <- RunTSNE(object = e16epithelium, dims = 1:10)
DimPlot(object = e16epithelium, reduction = "tsne", pt.size = 1.5,label = T,label.size = 6, cols = c("forestgreen", "navy", "red", "goldenrod2", "royalblue2", "blue",  "mediumpurple2", "darkorange3", "green3"))

### Identify E16 Unsupervised cluster markers
e16epithelium <- SetIdent(e16epithelium, value = "seurat_clusters")
e16epi.markers <- FindAllMarkers(object = e16epithelium, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.25)
e16epi.markers <- e16epi.markers[e16epi.markers$p_val_adj<0.05,]
write.table(e16epi.markers, file = "E16 Unsupervised Epithelium markers.txt", sep = "\t", row.names = F)
e16epi.markers.top5 <- e16epi.markers %>% group_by(cluster) %>% top_n(-5,p_val_adj)
DotPlot(object = e16epithelium, features = c(rev(unique(e16epi.markers.top5$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(angle = 90, size = 14))

### Generate tSNE plots for Figure
plot1 <- FeaturePlot(e16epithelium, features = c("Krt5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(e16epithelium, features = c("Aqp5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot3 <- FeaturePlot(e16epithelium, features = c("Krt14"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot4 <- FeaturePlot(e16epithelium, features = c("Acta2"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

plot1 <- FeaturePlot(e16epithelium, features = c("Ung"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(e16epithelium, features = c("Aurka"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(e16epithelium, features = c("Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)

plot1 <- FeaturePlot(e16epithelium, features = c("Igfbp2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(e16epithelium, features = c("Smgc"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot3 <- FeaturePlot(e16epithelium, features = c("Anxa1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot4 <- FeaturePlot(e16epithelium, features = c("Smoc2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot5 <- FeaturePlot(e16epithelium, features = c("Bpifa2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot6 <- FeaturePlot(e16epithelium, features = c("Elf5"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot7 <- FeaturePlot(e16epithelium, features = c("Barx1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot8 <- FeaturePlot(e16epithelium, features = c("Tubb5"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7,plot8), ncol = 1)

plot1 <- FeaturePlot(e16epithelium, features = c("Col4a2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(e16epithelium, features = c("Cnn1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot3 <- FeaturePlot(e16epithelium, features = c("Anxa1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))

## Annotate E16 epithelial clusters
e16epithelium <- RenameIdents(e16epithelium, '0' = "Myoepithelial", '1' = "End bud",'2'= "Krt19+ duct",'3'= "Basal duct",'4'= "End bud",'5'= "End bud",'6'= "Mesenchyme", '7' = "Basal duct", '8'= "Myoepithelial")
Idents(object = e16epithelium) <- factor(x = Idents(object = e16epithelium), levels = c("End bud", "Krt19+ duct", "Basal duct", "Myoepithelial", "Mesenchyme"))
e16epithelium[["Epi.CellType"]] <- Idents(object = e16epithelium)
e16epithelium <- SetIdent(e16epithelium, value = "Epi.CellType")
DimPlot(object = e16epithelium, reduction = "tsne", pt.size = 1.5,label = T,label.size = 8, cols = c("skyblue3", "red", "forestgreen", "hotpink", "goldenrod2")) + theme(axis.text = element_text(size = 14, colour = "black"))

### Determine differentially expressed genes in E16 epitheial clusters
e16.celltype.markers <- FindAllMarkers(e16epithelium,logfc.threshold = 0.25,only.pos = T)
e16.celltype.markers <- e16.celltype.markers[e16.celltype.markers$p_val_adj<0.05,]
write.table(e16.celltype.markers, file = "E16 cell type (SEURAT v3).txt", sep = "\t", row.names = F)
e16.celltype.markers.top10 <- e16.celltype.markers %>% group_by(cluster) %>% top_n(-10,p_val_adj)
DotPlot(e16epithelium, features = c(rev(unique(e16.celltype.markers.top10$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(size = 14))

## Belinda's markers
DotPlot(e16epithelium, group.by = "Epi.CellType", features = rev(c("Sema3b", "Ccl20", "Mgp", "Col9a1", "Id4", "Dcpp1", "Tinagl1", "Myl9", "Etv5", "Arhgdib", "Foxq1", "Ppp1r1b", "Nupr1", "Krt7", "Cyp2f2", "Tacstd2", "Wfdc18", "Runx1", "Krt19", "Wnt5b")), cols = c("Spectral")) + theme(axis.text.x =  element_text(angle = 90))

saveRDS(e16epithelium, file = "E16 Epithelium annotated (SEURAT v3).rds")


#######################################################################
##########################  P1  #######################################
#############  Re-cluster epithelium with SEURAT 

p1epithelium <- NormalizeData(object = p1epithelium, normalization.method = "LogNormalize", scale.factor = 10000)
p1epithelium <- FindVariableFeatures(object = p1epithelium, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = p1epithelium), 10)
plot1 <- VariableFeaturePlot(object = p1epithelium)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(x = p1epithelium)
p1epithelium <- ScaleData(object = p1epithelium, features = all.genes)
p1epithelium <- RunPCA(object = p1epithelium, features = VariableFeatures(object = p1epithelium))
ElbowPlot(object = p1epithelium)
p1epithelium <- FindNeighbors(object = p1epithelium, dims = 1:10, force.recalc = T)
p1epithelium <- FindClusters(object = p1epithelium, resolution = 0.8)
p1epithelium <- RunTSNE(object = p1epithelium, dims = 1:10)

DimPlot(object = p1epithelium, reduction = "tsne", pt.size = 1.5,label = T,label.size = 8, cols = c("plum", "purple", "turquoise3", "dodgerblue",  "tomato",  "goldenrod1", "hotpink",  "mediumpurple1", "tomato3",  "blue", "midnightblue",  "forestgreen","tan3", "darkred")) + theme(axis.text = element_text(size = 18, colour = "black"), axis.title =  element_text(size = 18, colour = "black")) + NoLegend()

### Identify P1 Unsupervised cluster markers
p1epithelium <- SetIdent(p1epithelium, value = "seurat_clusters")
p1epi.markers <- FindAllMarkers(object = p1epithelium, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.25)
p1epi.markers <- p1epi.markers[p1epi.markers$p_val_adj<0.05,]
write.table(p1epi.markers, file = "P1 Unsupervised Epithelium markers.txt", sep = "\t", row.names = F)
p1epi.markers.top5 <- p1epi.markers %>% group_by(cluster) %>% top_n(-5,p_val_adj)
p1epithelium <- SetIdent(p1epithelium, value = "seurat_clusters")
DotPlot(object = p1epithelium, features = c(rev(unique(p1epi.markers.top5$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(angle = 90, size = 14))

## Violin and Tsne plots for ID cell populations
FeaturePlot(p1epithelium, features = c("Krt14", "Krt19", "Krt5", "Smgc","Bpifa2", "Bhlha15", "Aqp5", "Cnn1", "Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  
### Generate tSNE plots for Figure
plot1 <- FeaturePlot(p1epithelium, features = c("Ssr4"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(p1epithelium, features = c("Smgc"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot3 <- FeaturePlot(p1epithelium, features = c("Bpifa2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot4 <- FeaturePlot(p1epithelium, features = c("Wfdc18"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot5 <- FeaturePlot(p1epithelium, features = c("Fbln2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot6 <- FeaturePlot(p1epithelium, features = c("Fxyd2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot7 <- FeaturePlot(p1epithelium, features = c("Lpo"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot8 <- FeaturePlot(p1epithelium, features = c("Cnn1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot9 <- FeaturePlot(p1epithelium, features = c("Cenpf"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7,plot8,plot9), ncol = 1)

plot1 <- FeaturePlot(p1epithelium, features = c("Ung"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(p1epithelium, features = c("Aurka"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(p1epithelium, features = c("Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)

plot1 <- FeaturePlot(p1epithelium, features = c("Aqp5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(p1epithelium, features = c("Bhlha15"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(p1epithelium, features = c("Smgc"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)

plot1 <- FeaturePlot(p1epithelium, features = c("Acta2"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(p1epithelium, features = c("Cnn1"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(p1epithelium, features = c("Myh11"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)

plot1 <- FeaturePlot(p1epithelium, features = c("Gfra3"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(p1epithelium, features = c("Kit"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(p1epithelium, features = c("Nkd2"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot2 <- FeaturePlot(p1epithelium, features = c("Esp18"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(p1epithelium, features = c("Prol1"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)

### P1 Duct starts to express markers of intercalated and striated ducts
tsne.list <- list()
features <- c("Krt19", "Klk1", "Ascl3", "Kit", "Nkd2", "Mki67", "Smgc", "Bpifa2")
for (i in 1:8){
  tsne.list[[i]] <-FeaturePlot(p1epithelium, features = features[i], pt.size = 1, min.cutoff = "q10", max.cutoff = "q90") + theme(axis.title = element_blank())
}
CombinePlots(plots = tsne.list, ncol = 1)

## Annotate P1 epithelial clusters
p1epithelium <- RenameIdents(p1epithelium, '0' = "Smgc+ Proacinar", '1' = "Smgc+ Proacinar",'2'= "Bpifa2+ Proacinar",'3'= "Bpifa2+ Proacinar",
                             '4'= "Duct cluster 4",'5'= "Basal duct",'6'= "Smgc+ Proacinar", '7' = "Smgc+ Proacinar", '8'= "Krt19+ duct",
                             '9'= "Bpifa2+ Proacinar",'10'= "Bpifa2+ Smgc+ Proacinar",'11'= "Myoepithelial",'12'= "Basal duct",'13'= "Krt19+ duct")

Idents(object = p1epithelium) <- factor(x = Idents(object = p1epithelium), levels = c("Smgc+ Proacinar", "Bpifa2+ Proacinar", "Bpifa2+ Smgc+ Proacinar", "Krt19+ duct", "Duct cluster 4", "Basal duct", "Myoepithelial"))
p1epithelium[["Epi.CellType"]] <- Idents(object = p1epithelium)
p1epithelium <- SetIdent(p1epithelium, value = "Epi.CellType")
DimPlot(object = p1epithelium, reduction = "tsne", pt.size = 1.5,label = T,label.size = 8, cols = c("purple","navy", "slateblue2", "darkred","firebrick2", "goldenrod2", "forestgreen")) + theme(axis.text = element_text(size = 14, colour = "black"))

### Determine differentially expressed genes in P1 epitheial clusters

p1.celltype.markers <- FindAllMarkers(p1epithelium,logfc.threshold = 0.25,only.pos = T)
p1.celltype.markers <- p1.celltype.markers[p1.celltype.markers$p_val_adj<0.05,]
write.table(p1.celltype.markers, file = "P1 cell type markers (SEURAT v3).txt", sep = "\t", row.names = F)
p1.celltype.markers.top10 <- p1.celltype.markers %>% group_by(cluster) %>% top_n(-10,p_val_adj)
DotPlot(p1epithelium, features = c(rev(unique(c("Epcam",p1.celltype.markers.top10$gene)))), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(size = 14))

saveRDS(p1epithelium, file = "P1 Epithelium annotated (SEURAT v3).rds")


############################################################################
########################## Adult epithelium  ###############################
#############  Re-cluster epithelium with SEURAT 
adultepithelium <- NormalizeData(object = adultepithelium, normalization.method = "LogNormalize", scale.factor = 10000)
adultepithelium <- FindVariableFeatures(object = adultepithelium, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x = adultepithelium)
adultepithelium <- ScaleData(object = adultepithelium, features = all.genes)
adultepithelium <- RunPCA(object = adultepithelium, features = VariableFeatures(object = adultepithelium))
ElbowPlot(object = adultepithelium)
adultepithelium <- FindNeighbors(object = adultepithelium, dims = 1:12, force.recalc = T)
adultepithelium <- FindClusters(object = adultepithelium, resolution = 0.6)
adultepithelium <- RunTSNE(object = adultepithelium, dims = 1:12)

DimPlot(object = adultepithelium, reduction = "tsne", pt.size = 1,label = T,label.size = 6, group.by = "seurat_clusters") + scale_color_manual(values = c("pink", "dodgerblue", "blue", "plum", "turquoise", "steelblue", "slateblue2", "purple", "goldenrod2", "orange", "firebrick2", "lightsteelblue", "forestgreen", "hotpink", "darkred", "navy"))
DimPlot(object = adultepithelium, reduction = "tsne", pt.size = 1,label = T,label.size = 6) 

  
### Identify Adult Unsupervised cluster markers
adultepithelium <- SetIdent(adultepithelium, value = "seurat_clusters")
adultepithelium.markers <- FindAllMarkers(object = adultepithelium, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.25)
adultepithelium.markers <- adultepithelium.markers[adultepithelium.markers$p_val_adj<0.05,]
write.table(adultepithelium.markers, file = "Adult Unsupervised Epithelium markers.txt", sep = "\t", row.names = F)

adultepi.markers.top5 <- adultepithelium.markers %>% group_by(cluster) %>% top_n(-5, p_val_adj)
DotPlot(object = adultepithelium, features = c(rev(unique(adultepi.markers.top5$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(angle = 90, size = 14))

## Violin and Tsne plots for ID cell populations
FeaturePlot(adultepithelium, features = c("Krt14", "Krt19", "Krt5", "Smgc","Bpifa2", "Bhlha15", "Aqp5", "Cnn1", "Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  
### Generate tSNE plots for Figure
plot1 <- FeaturePlot(adultepithelium, features = c("Smgc"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(adultepithelium, features = c("Lpo"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot3 <- FeaturePlot(adultepithelium, features = c("Fxyd2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot4 <- FeaturePlot(adultepithelium, features = c("Bpifa2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot5 <- FeaturePlot(adultepithelium, features = c("Nkd2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot6 <- FeaturePlot(adultepithelium, features = c("Ldhb"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot7 <- FeaturePlot(adultepithelium, features = c("Smoc2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot8 <- FeaturePlot(adultepithelium, features = c("Prol1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot9 <- FeaturePlot(adultepithelium, features = c("Cnn1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot10 <- FeaturePlot(adultepithelium, features = c("Ascl3"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot11 <- FeaturePlot(adultepithelium, features = c("S100a6"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot12 <- FeaturePlot(adultepithelium, features = c("Muc19"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))

CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5, plot6), ncol = 1)
CombinePlots(plots = list(plot7, plot8, plot9, plot10, plot11, plot12), ncol = 1)

plot1 <- FeaturePlot(adultepithelium, features = c("Ung"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(adultepithelium, features = c("Aurka"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(adultepithelium, features = c("Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)

plot1 <- FeaturePlot(adultepithelium, features = c("Aqp5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(adultepithelium, features = c("Bhlha15"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
CombinePlots(plots = list(plot1, plot2), ncol = 2)

DimPlot(adultepithelium, group.by = "CellType", label = T, label.size = 7)
## Annotate adult epithelial clusters
adultepithelium <- RenameIdents(adultepithelium, '0' = "Smgc+", '1' = "Acinar",'2'= "Acinar",'3'= "Smgc+",
                             '4'= "Striated duct",'5'= "Acinar",'6'= "Bpifa2+ Acinar", '7' = "Intercalated duct", '8'= "GCT",
                             '9'= "GCT",'10'= "Basal duct",'11'= "Acinar",'12'= "Myoepithelial",'13'= "Ascl3",'14'= "Basal duct", '15'="Bpifa2+ Acinar")

Idents(object = adultepithelium) <- factor(x = Idents(object = adultepithelium), levels = c("Bpifa2+ Acinar", "Acinar", "Smgc+", "Intercalated duct", "Basal duct","Striated duct", "Ascl3", "GCT", "Myoepithelial"))
adultepithelium[["Epi.CellType"]] <- Idents(object = adultepithelium)
adultepithelium <- SetIdent(adultepithelium, value = "Epi.CellType")
DimPlot(object = adultepithelium, reduction = "tsne", pt.size = 1.5,label = T,label.size = 6,repel = T, cols = c("slateblue2", "navy", "pink2", "purple", "firebrick2", "turquoise", "hotpink", "goldenrod2", "forestgreen")) + theme(axis.text = element_text(size = 14, colour = "black"))

### Determine differentially expressed genes in Adult epitheial clusters
adult.celltype.markers <- FindAllMarkers(adultepithelium,logfc.threshold = 0.25,only.pos = T)
adult.celltype.markers <- adult.celltype.markers[adult.celltype.markers$p_val_adj<0.05,]
write.table(adult.celltype.markers, file = "adult cell type markers (SEURAT v3).txt", sep = "\t", row.names = F)
adult.celltype.markers.toadult0 <- adult.celltype.markers %>% group_by(cluster) %>% top_n(-10,p_val_adj)

adultepithelium <- SetIdent(adultepithelium, value = "Epi.CellType")
DotPlot(adultepithelium, features = c(rev(unique(adult.celltype.markers.toadult0$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(size = 14))

plot1 <- VlnPlot(adultepithelium, idents = c("Acinar", "Bpifa2+ Acinar", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Bpifa2"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot2 <-VlnPlot(adultepithelium, idents = c("Acinar", "Bpifa2+ Acinar", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Bhlha15"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot3 <-VlnPlot(adultepithelium, idents = c("Acinar", "Bpifa2+ Acinar", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Smgc"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot4 <-VlnPlot(adultepithelium, idents = c("Acinar", "Bpifa2+ Acinar", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Gfra3"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot5 <-VlnPlot(adultepithelium, idents = c("Acinar", "Bpifa2+ Acinar", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Kit"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot6 <-VlnPlot(adultepithelium, idents = c("Acinar", "Bpifa2+ Acinar", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Krt7"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot7 <-VlnPlot(adultepithelium, idents = c("Acinar", "Bpifa2+ Acinar", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Prol1"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot8 <-VlnPlot(adultepithelium, idents = c("Acinar", "Bpifa2+ Acinar", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Aqp5"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8), ncol = 1)

saveRDS(adultepithelium, file = "Adult (C3H ICR) Epithelium annotated (SEURAT v3).rds")

############################################################################
########################## Adult epithelium  ###############################
########################## Integrated file (Seurat v2)  ####################
#############  Re-cluster epithelium with SEURAT 

DimPlot(adult.integrated, label = T, pt.size = 1)
adult.integrated.epithelium <- subset(adult.integrated,idents = c("Smgc+","Acinar","Bpifa2+", "Striated duct", "Basal duct", "Myoepithelial", "Ascl3+", "Intercalated duct", "GCT"))
DimPlot(adult.integrated.epithelium, label = T, pt.size = 1)

adult.integrated.epithelium <- NormalizeData(object = adult.integrated.epithelium, normalization.method = "LogNormalize", scale.factor = 10000)
adult.integrated.epithelium <- FindVariableFeatures(object = adult.integrated.epithelium, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x = adult.integrated.epithelium)
adult.integrated.epithelium <- ScaleData(object = adult.integrated.epithelium, features = all.genes)
adult.integrated.epithelium <- RunPCA(object = adult.integrated.epithelium, features = VariableFeatures(object = adult.integrated.epithelium))
ElbowPlot(object = adult.integrated.epithelium)
adult.integrated.epithelium <- FindNeighbors(object = adult.integrated.epithelium, dims = 1:10, force.recalc = T)
adult.integrated.epithelium <- FindClusters(object = adult.integrated.epithelium, resolution = 0.8)
adult.integrated.epithelium <- RunTSNE(object = adult.integrated.epithelium, dims = 1:10)

DimPlot(object = adult.integrated.epithelium, reduction = "tsne", pt.size = 1.5,label = T,label.size = 8, group.by = "seurat_clusters") + scale_color_manual(values = c("pink", "dodgerblue", "blue", "plum", "turquoise", "steelblue", "slateblue2", "purple", "goldenrod2", "orange", "firebrick2", "lightsteelblue", "forestgreen", "hotpink", "darkred", "navy"))
DimPlot(object = adult.integrated.epithelium, reduction = "tsne", pt.size = 1,label = T,label.size = 6) 

adult.integrated.epithelium <- SetIdent(adult.integrated.epithelium, value = "seurat_clusters")
adult.integrated.epithelium.markers <- FindAllMarkers(object = adult.integrated.epithelium, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.25)
adult.integrated.epithelium.markers <- adult.integrated.epithelium.markers[adult.integrated.epithelium.markers$p_val_adj<0.05,]
write.table(adult.integrated.epithelium.markers, file = "Adult (Integrated) Unsupervised Epithelium markers.txt", sep = "\t", row.names = F)

adult.integrated.unsup.markers <- adult.integrated.epithelium.markers %>% group_by(cluster) %>% top_n(-5, p_val_adj)
DotPlot(object = adult.integrated.epithelium, features = c(rev(unique(adult.integrated.unsup.markers$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(angle = 90, size = 14))

## Violin and Tsne plots for ID cell populations
FeaturePlot(adult.integrated.epithelium, features = c("Krt14", "Krt19", "Krt5", "Smgc","Bpifa2", "Bhlha15", "Aqp5", "Cnn1", "Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  
### Generate tSNE plots for Figure
plot1 <- FeaturePlot(adult.integrated.epithelium, features = c("Smgc"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(adult.integrated.epithelium, features = c("Lpo"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot3 <- FeaturePlot(adult.integrated.epithelium, features = c("Prol1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot4 <- FeaturePlot(adult.integrated.epithelium, features = c("Fxyd2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot5 <- FeaturePlot(adult.integrated.epithelium, features = c("Wfdc18"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot6 <- FeaturePlot(adult.integrated.epithelium, features = c("Krt14"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot9 <- FeaturePlot(adult.integrated.epithelium, features = c("Cnn1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot10 <- FeaturePlot(adult.integrated.epithelium, features = c("Ascl3"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot11 <- FeaturePlot(adult.integrated.epithelium, features = c("Bpifa2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))

CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5, plot6, plot9,plot10, plot11 ), ncol = 1)

FeaturePlot(adult.integrated.epithelium, features = c("Gfra3", "Smgc", "Kit", "Prol1"), pt.size = 0.5,min.cutoff = "q20", max.cutoff = "q90", cols = c("ghostwhite", "blue"), ncol = 1) + NoLegend()  + theme(title = element_blank(), axis.title = element_blank())
plot2 <- FeaturePlot(adult.integrated.epithelium, features = c("Smgc"), pt.size = 0.5,min.cutoff = "q20", max.cutoff = "q90", cols = c("ghostwhite", "blue")) + NoLegend()  + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot3 <- FeaturePlot(adult.integrated.epithelium, features = c("Kit"), pt.size = 0.5,min.cutoff = "q20", max.cutoff = "q90", cols = c("ghostwhite", "blue")) + NoLegend()  + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot4 <- FeaturePlot(adult.integrated.epithelium, features = c("Prol1"), pt.size = 0.5,min.cutoff = "q20", max.cutoff = "q90", cols = c("ghostwhite", "blue")) + NoLegend()  + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)


DimPlot(adult.integrated.epithelium, group.by = "seurat_clusters", label = T, label.size = 8, repel = F, pt.size = 1.5, cols = c("pink2", "dodgerblue", "blue", "turquoise", "purple2", "firebrick", "goldenrod2", "steelblue2", "forestgreen", "hotpink", "steelblue")) + theme(axis.text = element_text(size = 18, colour = "black"), axis.title =  element_text(size = 18, colour = "black")) + NoLegend()

## Annotate adult epithelial clusters
adult.integrated.epithelium <- RenameIdents(adult.integrated.epithelium, '0' = "Smgc+", '1' = "Acinar",'2'= "Acinar",
                                '3'= "Striated duct",'4' = "Intercalated duct", '5'="Basal duct", '6'= "GCT", '7'= "Bpifa2+", 
                                '8'= "Myoepithelial", '9'= "Ascl3+",'10'= "Bpifa2+")

Idents(object = adult.integrated.epithelium) <- factor(x = Idents(object = adult.integrated.epithelium), levels = c("Bpifa2+", "Acinar", "Smgc+", "Intercalated duct", "Basal duct","Striated duct", "Ascl3+", "GCT", "Myoepithelial"))
adult.integrated.epithelium[["Epi.CellType"]] <- Idents(object = adult.integrated.epithelium)
adult.integrated.epithelium <- SetIdent(adult.integrated.epithelium, value = "Epi.CellType")
DimPlot(object = adult.integrated.epithelium, reduction = "tsne", pt.size = 1.5,label = F,label.size = 6,repel = T, cols = c("slateblue2", "navy", "pink2", "purple", "firebrick2", "turquoise", "hotpink", "goldenrod2", "forestgreen")) + theme(axis.text = element_text(size = 14, colour = "black"))

### Determine differentially expressed genes in Adult epitheial clusters
adult.INT.celltype.markers <- FindAllMarkers(adult.integrated.epithelium,logfc.threshold = 0.25,only.pos = T)
adult.INT.celltype.markers <- adult.INT.celltype.markers[adult.INT.celltype.markers$p_val_adj<0.05,]
write.table(adult.INT.celltype.markers, file = "Adult (Integrated) cell type markers (SEURAT v3).txt", sep = "\t", row.names = F)
adult.INT.celltype.markers.top <- adult.INT.celltype.markers %>% group_by(cluster) %>% top_n(-10,p_val_adj)

adult.integrated.epithelium <- SetIdent(adult.integrated.epithelium, value = "Epi.CellType")

DotPlot(adult.integrated.epithelium, features = c(rev(unique(adult.INT.celltype.markers.top$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(size = 14))

plot1 <- VlnPlot(adult.integrated.epithelium, idents = c("Acinar", "Bpifa2+", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Aqp5"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot2 <-VlnPlot(adult.integrated.epithelium, idents = c("Acinar", "Bpifa2+", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Prol1"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot3 <-VlnPlot(adult.integrated.epithelium, idents = c("Acinar", "Bpifa2+", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Bhlha15"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot4 <-VlnPlot(adult.integrated.epithelium, idents = c("Acinar", "Bpifa2+", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Bpifa2"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot5 <-VlnPlot(adult.integrated.epithelium, idents = c("Acinar", "Bpifa2+", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Krt18"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot6 <-VlnPlot(adult.integrated.epithelium, idents = c("Acinar", "Bpifa2+", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Smgc"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot7 <-VlnPlot(adult.integrated.epithelium, idents = c("Acinar", "Bpifa2+", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Kit"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot8 <-VlnPlot(adult.integrated.epithelium, idents = c("Acinar", "Bpifa2+", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Gfra3"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot9 <-VlnPlot(adult.integrated.epithelium, idents = c("Acinar", "Bpifa2+", "Smgc+", "Intercalated duct"), cols = c("slateblue2", "blue", "pink2", "purple"), pt.size = 0,features = c("Krt7"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9), ncol = 1)

adult.integrated.epithelium <- SetIdent(adult.integrated.epithelium, value = "seurat_clusters")
plot1 <- VlnPlot(adult.integrated.epithelium, idents = c(1,2, 7,10, 0, 4), cols = c("pink2", "dodgerblue", "blue", "purple2", "steelblue2", "steelblue"), pt.size = 0,features = c("Aqp5"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot2 <-VlnPlot(adult.integrated.epithelium, idents = c(1,2, 7,10, 0, 4), cols = c("pink2", "dodgerblue", "blue", "purple2", "steelblue2", "steelblue"), pt.size = 0,features = c("Prol1"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot3 <-VlnPlot(adult.integrated.epithelium, idents = c(1,2, 7,10, 0, 4), cols = c("pink2", "dodgerblue", "blue", "purple2", "steelblue2", "steelblue"), pt.size = 0,features = c("Bhlha15"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot4 <-VlnPlot(adult.integrated.epithelium, idents = c(1,2, 7,10, 0, 4), cols = c("pink2", "dodgerblue", "blue", "purple2", "steelblue2", "steelblue"), pt.size = 0,features = c("Bpifa2"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot5 <-VlnPlot(adult.integrated.epithelium, idents = c(1,2, 7,10, 0, 4), cols = c("pink2", "dodgerblue", "blue", "purple2", "steelblue2", "steelblue"), pt.size = 0,features = c("Krt18"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot6 <-VlnPlot(adult.integrated.epithelium, idents = c(1,2, 7,10, 0, 4), cols = c("pink2", "dodgerblue", "blue", "purple2", "steelblue2", "steelblue"), pt.size = 0,features = c("Smgc"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot7 <-VlnPlot(adult.integrated.epithelium, idents = c(1,2, 7,10, 0, 4), cols = c("pink2", "dodgerblue", "blue", "purple2", "steelblue2", "steelblue"), pt.size = 0,features = c("Kit"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot8 <-VlnPlot(adult.integrated.epithelium, idents = c(1,2, 7,10, 0, 4), cols = c("pink2", "dodgerblue", "blue", "purple2", "steelblue2", "steelblue"), pt.size = 0,features = c("Gfra3"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
plot9 <-VlnPlot(adult.integrated.epithelium, idents = c(1,2, 7,10, 0, 4), cols = c("pink2", "dodgerblue", "blue", "purple2", "steelblue2", "steelblue"), pt.size = 0,features = c("Krt7"), ncol = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 12)) + NoLegend()
CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9), ncol = 1)


saveRDS(adult.integrated.epithelium, file = "Adult (Integrated C3H ICR) Epithelium annotated (SEURAT v3).rds")

plot1 <- FeaturePlot(adult.integrated.epithelium, features = c("Aqp5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(adult.integrated.epithelium, features = c("Bhlha15"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(adult.integrated.epithelium, features = c("Smgc"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)

plot1 <- FeaturePlot(adult.integrated.epithelium, features = c("Ung"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(adult.integrated.epithelium, features = c("Aurka"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(adult.integrated.epithelium, features = c("Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)

plot1 <- FeaturePlot(adult.integrated.epithelium, features = c("Col1a1"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
plot2 <- FeaturePlot(adult.integrated.epithelium, features = c("Col1a2"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
plot3 <- FeaturePlot(adult.integrated.epithelium, features = c("Barx1"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)

FeaturePlot(adult.integrated.epithelium, features = c("Gfra3"),cols = c("ghostwhite", "blue"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend()
FeaturePlot(adult.integrated.epithelium, features = c("Krt14"),cols = c("ghostwhite", "blue"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend()
FeaturePlot(adult.integrated.epithelium, features = c("Kit"),cols = c("ghostwhite", "blue"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend()
FeaturePlot(adult.integrated.epithelium, features = c("Gstt1"),cols = c("ghostwhite", "blue"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend()
FeaturePlot(adult.integrated.epithelium, features = c("Smgc"),cols = c("ghostwhite", "blue"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend()
FeaturePlot(adult.integrated.epithelium, features = c("Ret"),cols = c("ghostwhite", "blue"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend()
FeaturePlot(adult.integrated.epithelium, features = c("Prol1"),cols = c("ghostwhite", "blue"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend()

FeaturePlot(adult.integrated.epithelium, features = c("Gfra3"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt"))
FeaturePlot(adult.integrated.epithelium, features = c("Kit"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
FeaturePlot(adult.integrated.epithelium, features = c("Nkd2"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
FeaturePlot(adult.integrated.epithelium, features = c("Esp18"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
FeaturePlot(adult.integrated.epithelium, features = c("Prol1"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
FeaturePlot(adult.integrated.epithelium, features = c("Gstt1"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
FeaturePlot(adult.integrated.epithelium, features = c("Krt14"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )

#### Formatted plots for Figures #######

pdf(file = "Fig2-3 - E12 tSNEs.pdf", width = 5,height = 4)
DimPlot(object = e12epithelium, reduction = "tsne", pt.size = 1.5,label = F,group.by = "seurat_clusters", cols = c("skyblue1", "red", "darkred", "skyblue3", "royalblue2", "mediumpurple1", "hotpink", "gray")) + theme(axis.text = element_text(size = 18, colour = "black"), axis.title =  element_text(size = 18, colour = "black")) 
DimPlot(object = e14epithelium, group.by = "seurat_clusters", reduction = "tsne", pt.size = 1,label = F, cols = c("blue", "skyblue3", "royalblue2", "turquoise1", "mediumpurple2", "orange", "darkred")) + theme(axis.text = element_text(size = 18, colour = "black"), axis.title =  element_text(size = 18, colour = "black")) +NoLegend()
DimPlot(object = e16epithelium,group.by = "seurat_clusters", reduction = "tsne", pt.size = 1.2,label = F, cols = c("forestgreen", "navy", "red", "goldenrod2", "royalblue2", "blue",  "mediumpurple2", "darkorange3", "green3")) + theme(axis.text = element_text(size = 18, colour = "black"), axis.title =  element_text(size = 18, colour = "black")) 
DimPlot(object = p1epithelium, reduction = "tsne",group.by = "seurat_clusters", pt.size = 1,label = F,label.size = 8, cols = c("plum", "purple", "turquoise3", "dodgerblue",  "tomato",  "goldenrod1", "hotpink",  "mediumpurple1", "tomato3",  "blue", "midnightblue",  "forestgreen","tan3", "darkred")) + theme(axis.text = element_text(size = 18, colour = "black"), axis.title =  element_text(size = 18, colour = "black")) 
DimPlot(adult.integrated.epithelium, group.by = "seurat_clusters", label = F, repel = F, pt.size = 1, cols = c("pink2", "dodgerblue", "blue", "turquoise", "purple2", "firebrick", "goldenrod2", "steelblue2", "forestgreen", "hotpink", "steelblue")) + theme(axis.text = element_text(size = 18, colour = "black"), axis.title =  element_text(size = 18, colour = "black")) 
dev.off()

FeaturePlot(adult.integrated.epithelium, features = c("Dcdc2a"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
FeaturePlot(p1epithelium, features = c("Gstt1"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )
FeaturePlot(p1epithelium, features = c("Dcdc2a"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(1,1,0,1),units = "pt") )


e12epithelium <- SetIdent(e12epithelium, value = "seurat_clusters")
target <- c(0,3,4,1,2,5,6)
Idents(e12epithelium) <- factor(Idents(e12epithelium), levels = target)
require(gdata)
e12epi.markers.top5$cluster <- factor(e12epi.markers.top5$cluster, levels = target)
e12epi.markers.top5 <- e12epi.markers.top5[order(e12epi.markers.top5$cluster), ]
DotPlot(e12epithelium, features = c(rev(unique(e12epi.markers.top5$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(size = 14))

e14epithelium <- SetIdent(e14epithelium, value = "seurat_clusters")
DimPlot(e14epithelium, label = T)
target <- c(0,1,2,3,6,5,4)
Idents(e14epithelium) <- factor(Idents(e14epithelium), levels = target)
require(gdata)
e14epi.markers.top5$cluster <- factor(e14epi.markers.top5$cluster, levels = target)
e14epi.markers.top5 <- e14epi.markers.top5[order(e14epi.markers.top5$cluster), ]
DotPlot(e14epithelium, features = c(rev(unique(e14epi.markers.top5$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(size = 14))

e16epithelium <- SetIdent(e16epithelium, value = "seurat_clusters")
target <- c(1,4,5,2,3,7,0,8,6)
Idents(e16epithelium) <- factor(Idents(e16epithelium), levels = target)
require(gdata)
e16epi.markers.top5$cluster <- factor(e16epi.markers.top5$cluster, levels = target)
e16epi.markers.top5 <- e16epi.markers.top5[order(e16epi.markers.top5$cluster), ]
DotPlot(e16epithelium, features = c(rev(unique(e16epi.markers.top5$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(size = 14))

p1epithelium <- SetIdent(p1epithelium, value = "seurat_clusters")
target <- c(2,3,9,10,0,1,6,7,4,13,5,8,12,11)
Idents(p1epithelium) <- factor(Idents(p1epithelium), levels = target)
require(gdata)
p1epi.markers.top5$cluster <- factor(p1epi.markers.top5$cluster, levels = target)
p1epi.markers.top5 <- p1epi.markers.top5[order(p1epi.markers.top5$cluster), ]
DotPlot(p1epithelium, features = c(rev(unique(p1epi.markers.top5$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(size = 14))

adult.integrated.epithelium <- SetIdent(adult.integrated.epithelium, value = "seurat_clusters")
target <- c(1,2,10,7,0,4,3,6,9,5,8)
Idents(adult.integrated.epithelium) <- factor(Idents(adult.integrated.epithelium), levels = target)
require(gdata)
adult.integrated.unsup.markers$cluster <- factor(adult.integrated.unsup.markers$cluster, levels = target)
adult.integrated.unsup.markers <- adult.integrated.unsup.markers[order(adult.integrated.unsup.markers$cluster), ]
DotPlot(adult.integrated.epithelium, features = c(rev(unique(adult.integrated.unsup.markers$gene)), "Epcam"), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), axis.text.y = element_text(size = 14))

AvgExp.all.cells <- AverageExpression(adult.integrated.epithelium, slot="scale.data")
write.csv(AvgExp.all.cells, file = "Avg expression all clusters in adult epithelium.csv")
gfra3.vs.smgc <- FindMarkers(adult.integrated.epithelium, ident.1 = "Smgc+", ident.2 = "Intercalated duct", only.pos = F)
gfra3.vs.smgc <- gfra3.vs.smgc[gfra3.vs.smgc$p_val_adj<0.05, ]
write.csv(gfra3.vs.smgc, file = "Gfra3 vs Smgc cells.csv")


##### GFRA3 and SMGC heatmaps ######
gfra3.features <- FindMarkers(adult.integrated.epithelium, ident.1 = "Intercalated duct", group.by = "Epi.CellType", only.pos = T)
gfra3.features$gene <- rownames(gfra3.features)
gfra3.features<-gfra3.features[gfra3.features$p_val_adj<0.05, ]
gfra3.features.top100 <- gfra3.features[1:100,]

adult.integrated.epithelium<- SetIdent(adult.integrated.epithelium, value = "Epi.CellType")
avgexpGfra3 <- AverageExpression(object = adult.integrated.epithelium,features = rownames(gfra3.vs.smgc$gene))

library(pheatmap)
library(viridis)

avgexpGfra3.RNA <- t(avgexpGfra3$RNA)
avgexpGfra3.int <- t(avgexpGfra3$integrated)
avgexpGfra3.RNA <- avgexpGfra3.RNA[c(c(1:3), c(5:9, 4)), ]
pheatmap(avgexpGfra3.int, fontsize=10, color = magma(6), breaks = c(-2:3))

gfra3.vs.smgc <- gfra3.vs.smgc[order(gfra3.vs.smgc$avg_logFC), ]
DoHeatmap(object = subset(adult.integrated.epithelium, idents = c("Smgc+", "Intercalated duct") ), features = rownames(gfra3.vs.smgc), label = F, disp.min = -1) + theme(axis.text = element_blank()) + scale_fill_viridis(option = "magma")

superheat(avgexpGfra3.int,
          scale = T,legend = T, pretty.order.cols = T,
          bottom.label.text.angle = 90, bottom.label.size = 0.5, bottom.label.text.alignment = "right", bottom.label.text.size = 4,
          left.label.text.size = 4, left.label.col = "white", bottom.label.col = "white", grid.hline.size = 0.2, grid.vline.size = 0.2)

adult.integrated.epithelium <- SetIdent(adult.integrated.epithelium, value = "Epi.CellType")
smgc.features <- FindMarkers(adult.integrated.epithelium, ident.1 = "Smgc+", group.by = "Epi.CellType", only.pos = T)
smgc.features$gene <- rownames(smgc.features)
smgc.features<-smgc.features[smgc.features$p_val_adj<0.05, ]
smgc.features.top100 <- smgc.features[1:100,]

avgexpSmgc <- AverageExpression(object = adult.integrated.epithelium,features = smgc.features$gene)

avgexpSmgc.RNA <- t(avgexpSmgc$RNA)
avgexpSmgc.int <- t(avgexpSmgc$integrated)

avgexpSmgc.RNA <- avgexpSmgc.RNA[c(c(3, 1:2), c(4:9)), ]

superheat(avgexpSmgc.RNA,
          scale = T,legend = T, pretty.order.cols = F,
          bottom.label.text.angle = 90, bottom.label.size = 0.5, bottom.label.text.alignment = "right", bottom.label.text.size = 4,
          left.label.text.size = 4, left.label.col = "white", bottom.label.col = "white", grid.hline.size = 0.2, grid.vline.size = 0.2)


superheat(avgexpSmgc.RNA,
          scale = T,legend = T, pretty.order.cols = F,
          bottom.label.text.angle = 90, bottom.label.size = 0.5, bottom.label.text.alignment = "right", bottom.label.text.size = 4,
          left.label.text.size = 4, left.label.col = "white", bottom.label.col = "white", grid.hline.size = 0.2, grid.vline.size = 0.2)


#### SEPARATE SMGC AND ID CELLS TO MERGE WITH TABULA MURIS FILES ####
DimPlot(adult.integrated.epithelium)
ID.cells <- subset(adult.integrated.epithelium, idents = c("Smgc+", "Intercalated duct"))
saveRDS(ID.cells, file = "ID cells subset.rds")
DimPlot(ID.cells)

test <- gfra3.features[gfra3.features$gene %in% smgc.features$gene, ]
test <- test[!(test$gene %in% rownames(gfra3.vs.smgc)), ]
