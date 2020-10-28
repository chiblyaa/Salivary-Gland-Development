# Script for analysis of embryonic SMG epithelium 
# associated with Figure 2 

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(superheat)
library(Hmisc)
library(clustree)
library(cowplot)

#use same color scheme as in previous script for consistency
colors <- c("#eb7341",
            "#8a308c",
            "#4ab539",
            "#ffb558",
            "#006f35",
            "#b65ea4",
            "#b8cd42",
            "#f3e3cb", 
            "#adf4d0",  
            "#565e00",
            "#47b5d5",
            "#ff9e96",
            "#01908c",
            "#7152f6",
            "#9a8e62",
            "#e10060",
            "#345a78",
            "#9c001b",
            "#eee799",
            "#e74134", 
            "#0074c9", 
            "#7f4a21", 
            "#c3c1ff",
            "#fbe54b")

cell.types <- c("Acinar", "Ascl3+ duct", "Basal duct", "Bpifa2+", "Bpifa2+ Proacinar", "End bud", "Endothelial",
                "Erythroid", "GCT", "Glial cells", "Intercalated duct", "Krt19+ duct", "Macrophages", "Mast cells",
                "Mesenchyme", "Mitotic cells", "Myoepithelial", "Nerves", "NK cells", "Smgc+", "Smgc+ Proacinar",
                "Smooth muscle", "Striated duct", "Stromal")

names(colors) <- cell.types

stage.colors <- c(
  "#c34d54",
  "#d7cd9e",
  "#019f84",
  "#92d8d4",
  "#0076af",
  "#6b0081"
)
names(stage.colors) <- c("E12", "E14", "E16", "P1", "P30", "Adult")

########################## Load previously annotated SEURAT files############################
#############################################################################################

e.epithelium <- readRDS("../Embryonic epithelium.rds")
DefaultAssay(e.epithelium) <- "RNA"
pdf(file = "Embryonic epithelium (subset).pdf", useDingbats = F, width = 5.5, height = 4)
DimPlot(object = e.epithelium, pt.size = 1,label = F,repel = F, label.size = 5, group.by = "CellType", cols = colors) + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank())
dev.off()

### Split object to re-cluster epithelium from individual developmental stages
e.split <- SplitObject(e.epithelium, split.by = "stage")
### standard worflow for normalization and scaling
for (i in 1:length(e.split)) {
  e.split[[i]] <- NormalizeData(object = e.split[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  e.split[[i]] <- FindVariableFeatures(object = e.split[[i]], selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x = e.split[[i]])
  e.split[[i]] <- ScaleData(object = e.split[[i]], features = all.genes, vars.to.regress = "percent.mt")
  e.split[[i]] <- RunPCA(object = e.split[[i]], features = VariableFeatures(object = e.split[[i]]))
}

# extract individual stages after normalization and scaling
e12epithelium <- e.split[[1]]
e14epithelium <- e.split[[2]]
e16epithelium <- e.split[[3]]

save(e12epithelium, e14epithelium, e16epithelium, file = "Embryonic epithelium re-scaled before reclustering (backup).RData")

######## ANALYSIS OF E12 epitheliumTHELIUM 
#Evaluate optimal clustering resolution using clustree package
DefaultAssay(e12epithelium) <- "RNA"
ElbowPlot(object = e12epithelium)
e12epithelium <- FindNeighbors(object = e12epithelium, dims = 1:10,force.recalc = T)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,0.9, 1, 1.2, 1.4)){
  e12epithelium <- FindClusters(e12epithelium, resolution = resolution)
}
pdf('E12_clustree.pdf', width = 15, height = 10)
clustree(e12epithelium,prefix = "RNA_snn_res.")
dev.off()
# optimal resolution at 0.5. At higher resolution, clusters from one branch have multiple parents which is abnormal.
# perform clustering
e12epithelium <- FindClusters(object = e12epithelium, resolution = 0.5)
e12epithelium <- RunUMAP(object = e12epithelium, dims = 1:10)

pdf(file = "E12 epithelium re-clustered (original label).pdf", useDingbats = F, width = 5.5, height = 4)
DimPlot(object = e12epithelium, pt.size = 1,label = T,repel = F, label.size = 5, group.by = "CellType", cols = colors) + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank())
dev.off()
pdf(file = "E12 epithelium re-clustered (new cluster IDs).pdf", useDingbats = F, width = 3, height = 3)
DimPlot(object = e12epithelium, pt.size = 1,label = T,repel = F, label.size = 5, group.by = "RNA_snn_res.0.5") + NoAxes() + NoLegend()
dev.off()

e12epithelium$unsupervised <- Idents(e12epithelium)

### Identify E12 markers for putative epithelial subpopulations
e12epithelium <- SetIdent(e12epithelium, value = "unsupervised")
e12epi.markers <- FindAllMarkers(e12epithelium, logfc.threshold = 0.25, only.pos = T)
write.csv(e12epi.markers, file = "E12 Unsupervised Epithelium markers.csv")

#Re-order identities for easier visualization of subclusters from similar populations
Idents(e12epithelium) <- factor(Idents(e12epithelium), levels = c(0,3,4,2,1,5,6))
DimPlot(e12epithelium)

target <- c(0,3,4,2,1,5,6)
e12epi.markers$cluster <- factor(e12epi.markers$cluster, levels = target)
e12epi.markers.top5 <- e12epi.markers %>% group_by(cluster) %>% top_n(-5,p_val_adj)

pdf(file = "E12 epithelium balloon plot.pdf", useDingbats = F, width = 9, height = 3)
DotPlot(object = e12epithelium, features = c("Epcam",unique(e12epi.markers.top5$gene)), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(angle = 90, size = 14, hjust = 0.5), axis.title = element_blank())
dev.off()

### Generate tSNE plots for Figure
library(cowplot)
plot1 <- FeaturePlot(e12epithelium, features = c("Krt5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(e12epithelium, features = c("Krt19"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16) )
plot3 <- FeaturePlot(e12epithelium, features = c("Krt14"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-0.4,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(e12epithelium, features = c("Sox10"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes()+ theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-0.4,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("E12 - Krt UMAPS Fig2.pdf", useDingbats = F, width = 5, height = 5)
plot_grid(plotlist = list(plot1,plot2,plot3,plot4), nrow = 2)
dev.off()

plot1 <- FeaturePlot(e12epithelium, features = c("Ung"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(e12epithelium, features = c("Aurka"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16) )
plot3 <- FeaturePlot(e12epithelium, features = c("Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-0.4,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("E12 - Proliferation UMAPS Fig2.pdf", useDingbats = F, width = 5, height = 5)
plot_grid(plotlist = list(plot1,plot2,plot3), nrow = 2)
dev.off()

plot1 <- FeaturePlot(e12epithelium, features = c("Snhg4"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16)) 
plot2 <- FeaturePlot(e12epithelium, features = c("Krt15"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot3 <- FeaturePlot(e12epithelium, features = c("Wnt6"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(e12epithelium, features = c("Ccne1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot5 <- FeaturePlot(e12epithelium, features = c("Ccnb1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot6 <- FeaturePlot(e12epithelium, features = c("Col1a2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot7 <- FeaturePlot(e12epithelium, features = c("Eif4a1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("E12 - Representative UMAPS Fig2.pdf", useDingbats = F, width = 2.5, height = 15)
plot_grid(plotlist = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7), ncol = 1)
dev.off()

### Annotate epithelial clusters
e12epithelium <- SetIdent(e12epithelium, value = "unsupervised")
e12epithelium <- RenameIdents(e12epithelium, '0' = "End bud", '1' = "Krt19+ duct", '2'= "Krt19+ duct", '3'= "End bud",'4'= "End bud",'5'= "Mesenchyme",'6'= "Undefined")
Idents(object = e12epithelium) <- factor(x = Idents(object = e12epithelium), levels = c("End bud", "Krt19+ duct", "Mesenchyme", "Undefined"))
e12epithelium[["celltype.fixed"]] <- Idents(object = e12epithelium)
e12epithelium <- SetIdent(e12epithelium, value = "celltype.fixed")

e12epithelium <- SetIdent(e12epithelium, value = "unsupervised")
e12epithelium <- RenameIdents(e12epithelium, '0' = "End bud", '1' = "Krt19+ duct", '2'= "Proliferative Krt19+ duct", '3'= "End bud",'4'= "Proliferative End bud",'5'= "Mesenchyme",'6'= "Undefined")
e12epithelium[["subpopulations"]] <- Idents(object = e12epithelium)
e12epithelium <- SetIdent(e12epithelium, value = "subpopulations")
DimPlot(e12epithelium)

pdf(file = "E12 epithelium re-clustered (Fixed CellType IDs).pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(object = e12epithelium, pt.size = 1.5,label = T,label.size = 5, repel = T, cols = colors, group.by = "celltype.fixed") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank()) + NoLegend()
dev.off()
pdf(file = "E12 epithelium re-clustered (Subpopulation IDs).pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(object = e12epithelium, pt.size = 1.5,label = T,label.size = 5, repel = T, group.by = "subpopulations") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank()) + NoLegend()
dev.off()

e12epithelium <- SetIdent(e12epithelium, value = "celltype.fixed")
e12cell.markers <- FindAllMarkers(e12epithelium, logfc.threshold = 0.25, only.pos = T)
e12cell.markers <- e12cell.markers[e14cell.markers$p_val_adj<0.05, ]
write.csv(e12cell.markers, file = "E12 Cell Type markers (Epithelium).csv")

saveRDS(e12epithelium, file = "E12 epithelium qcd_annotated.rds")

######## ANALYSIS OF E14 epithelium
#Evaluate optimal clustering resolution using clustree package
DefaultAssay(e14epithelium) <- "RNA"
ElbowPlot(object = e14epithelium)
e14epithelium <- FindNeighbors(object = e14epithelium, dims = 1:10,force.recalc = T)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,0.9, 1, 1.2, 1.4)){
  e14epithelium <- FindClusters(e14epithelium, resolution = resolution)
}
pdf('e14_clustree.pdf', width = 15, height = 10)
clustree(e14epithelium,prefix = "RNA_snn_res.")
dev.off()

# perform clustering
e14epithelium <- FindClusters(object = e14epithelium, resolution = 0.5)
e14epithelium <- RunUMAP(object = e14epithelium, dims = 1:10)
e14epithelium$unsupervised <- Idents(e14epithelium)

pdf(file = "e14 epithelium re-clustered (original label).pdf", useDingbats = F, width = 5.5, height = 4)
DimPlot(object = e14epithelium, pt.size = 1,label = T,repel = F, label.size = 5, group.by = "CellType", cols = colors) + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank())
dev.off()
pdf(file = "e14 epithelium re-clustered (new cluster IDs).pdf", useDingbats = F, width = 3, height = 3)
DimPlot(object = e14epithelium, pt.size = 1,label = T,repel = F, label.size = 5, group.by = "unsupervised") + NoLegend() +NoAxes()
dev.off()


### Identify E14 markers for putative epithelial subpopulations
e14epithelium <- SetIdent(e14epithelium, value = "unsupervised")
e14epi.markers <- FindAllMarkers(e14epithelium, logfc.threshold = 0.25, only.pos = T)
write.csv(e14epi.markers, file = "e14 Unsupervised Epithelium markers.csv")

#Re-order identities for easier visualization of subclusters from similar populations
Idents(e14epithelium) <- factor(Idents(e14epithelium), levels = c(0,1,2,3,5,6,4))
DimPlot(e14epithelium)

target <- c(0,1,2,3,5,6,4)
e14epi.markers$cluster <- factor(e14epi.markers$cluster, levels =  c(0,1,2,3,5,6,4))
e14epi.markers.top5 <- e14epi.markers %>% group_by(cluster) %>% top_n(-5,p_val_adj)
e14epi.markers.top5 <- e14epi.markers.top5[order(e14epi.markers.top5$cluster), ]

pdf(file = "e14 epithelium balloon plot.pdf", useDingbats = F, width = 9, height = 3)
DotPlot(object = e14epithelium, features = c("Epcam",unique(e14epi.markers.top5$gene)), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(angle = 90, size = 14, hjust = 0.5), axis.title = element_blank())
dev.off()


### Generate UMAP plots for Figure
plot1 <- FeaturePlot(e14epithelium, features = c("Krt5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(e14epithelium, features = c("Krt19"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16) )
plot3 <- FeaturePlot(e14epithelium, features = c("Krt14"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-0.4,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(e14epithelium, features = c("Sox10"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes()+ theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-0.4,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("e14 - Krt UMAPS Fig2.pdf", useDingbats = F, width = 5, height = 5)
plot_grid(plotlist = list(plot1,plot2,plot3,plot4), nrow = 2)
dev.off()

plot1 <- FeaturePlot(e14epithelium, features = c("Ung"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(e14epithelium, features = c("Aurka"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16) )
plot3 <- FeaturePlot(e14epithelium, features = c("Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-0.4,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("e14 - Proliferation UMAPS Fig2.pdf", useDingbats = F, width = 5, height = 5)
plot_grid(plotlist = list(plot1,plot2,plot3), nrow = 2)
dev.off()

plot1 <- FeaturePlot(e14epithelium, features = c("Actg1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.94, vjust = -8, size = 16)) 
plot2 <- FeaturePlot(e14epithelium, features = c("Ccne1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.94, vjust = -8, size = 16))
plot3 <- FeaturePlot(e14epithelium, features = c("Cenpf"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.94, vjust = -8, size = 16))
plot4 <- FeaturePlot(e14epithelium, features = c("Ccnd1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.94, vjust = -8, size = 16))
plot5 <- FeaturePlot(e14epithelium, features = c("Col1a2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.94, vjust = -8, size = 16))
plot6 <- FeaturePlot(e14epithelium, features = c("Foxq1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.94, vjust = -8, size = 16))
plot7 <- FeaturePlot(e14epithelium, features = c("Id3"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.94, vjust = -8, size = 16))
pdf("e14 - Representative UMAPS Fig2.pdf", useDingbats = F, width = 2.5, height = 15)
plot_grid(plotlist = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7), ncol = 1)
dev.off()

pdf("E14 - Cldn10 UMAP Fig2.pdf", useDingbats = F, width = 2.5, height = 2.5)
FeaturePlot(e14epithelium, features = c("Cldn10"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
dev.off()


## Annotate epithelial clusters
e14epithelium <- SetIdent(e14epithelium, value = "unsupervised")
e14epithelium <- RenameIdents(e14epithelium, '0' = "End bud", '1' = "End bud",'2'= "End bud",'3'= "End bud",'4'= "Mesenchyme",'5'= "Krt19+ duct",'6'= "Basal duct")
Idents(object = e14epithelium) <- factor(x = Idents(object = e14epithelium), levels = c("End bud", "Krt19+ duct", "Basal duct", "Mesenchyme"))
e14epithelium[["celltype.fixed"]] <- Idents(object = e14epithelium)
e14epithelium <- SetIdent(e14epithelium, value = "celltype.fixed")

e14epithelium <- SetIdent(e14epithelium, value = "unsupervised")
e14epithelium <- RenameIdents(e14epithelium, '0' = "End bud", '1' = "End bud", '2'= "Proliferative End bud", '3'= "Cldn10+ End bud",'4'= "Mesenchyme",'5'= "Krt19+ duct",'6'= "Basal duct")
e14epithelium[["subpopulations"]] <- Idents(object = e14epithelium)
e14epithelium <- SetIdent(e14epithelium, value = "subpopulations")
DimPlot(e14epithelium)

pdf(file = "e14 epithelium re-clustered (Fixed CellType IDs).pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(object = e14epithelium, pt.size = 1.5,label = T,label.size = 5, repel = T, cols = colors, group.by = "celltype.fixed") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank()) + NoLegend()
dev.off()
pdf(file = "e14 epithelium re-clustered (Subpopulation IDs).pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(object = e14epithelium, pt.size = 1.5,label = T,label.size = 5, repel = T, group.by = "subpopulations") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank()) + NoLegend()
dev.off()

pdf(file = "e14 epithelium re-clustered (mesenchymal markers).pdf", useDingbats = F, width = 5.5, height = 5)
FeaturePlot(e14epithelium, features = "Col1a1", min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(e14epithelium, features = "Col3a1", min.cutoff = "q10", max.cutoff = "q90")
dev.off()

e14epithelium <- SetIdent(e14epithelium, value = "celltype.fixed")
e14cell.markers <- FindAllMarkers(e14epithelium, logfc.threshold = 0.25, only.pos = T)
e14cell.markers <- e14cell.markers[e14cell.markers$p_val_adj<0.05, ]
write.csv(e14cell.markers, file = "E14 Cell Type markers (Epithelium).csv")

saveRDS(e14epithelium, file = "E14 epithelium qcd_annotated.rds")

######## ANALYSIS OF E16 epithelium
#Evaluate optimal clustering resolution using clustree package
DefaultAssay(e16epithelium) <- "RNA"
ElbowPlot(object = e16epithelium)
e16epithelium <- FindNeighbors(object = e16epithelium, dims = 1:12,force.recalc = T)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,0.9, 1, 1.2, 1.4)){
  e16epithelium <- FindClusters(e16epithelium, resolution = resolution)
}
pdf('e16_clustree.pdf', width = 15, height = 10)
clustree(e16epithelium,prefix = "RNA_snn_res.")
dev.off()
#optimal resolution for ID of cells was 0.3
#however, to identify subpopulations, a resolution of 0.8 was needed
# which identified subclusters of proliferative cells

# perform clustering
e16epithelium <- FindClusters(object = e16epithelium, resolution = 0.8)
e16epithelium <- RunUMAP(object = e16epithelium, dims = 1:12)
e16epithelium$unsupervised <- Idents(e16epithelium)

pdf(file = "e16 epithelium re-clustered (original label).pdf", useDingbats = F, width = 5.5, height = 4)
DimPlot(object = e16epithelium, pt.size = 1,label = T,repel = T, label.size = 5, group.by = "CellType", cols = colors) + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank())
dev.off()
pdf(file = "e16 epithelium re-clustered (new cluster IDs).pdf", useDingbats = F, width = 3, height = 3)
DimPlot(object = e16epithelium, pt.size = 1,label = T,repel = F, label.size = 5, group.by = "unsupervised") + NoLegend() + NoAxes()
dev.off()


### Identify E16 markers for putative epithelial subpopulations
e16epithelium <- SetIdent(e16epithelium, value = "unsupervised")
e16epi.markers <- FindAllMarkers(e16epithelium, logfc.threshold = 0.25, only.pos = T)
write.csv(e16epi.markers, file = "e16 Unsupervised Epithelium markers.csv")

#Re-order identities for easier visualization of subclusters from similar populations
target <- c(1,4,5,2,3,8,0,7,6)
Idents(e16epithelium) <- factor(Idents(e16epithelium), levels = target)
DimPlot(e16epithelium)

e16epi.markers$cluster <- factor(e16epi.markers$cluster, levels =  target)
e16epi.markers.top5 <- e16epi.markers %>% group_by(cluster) %>% top_n(-5,p_val_adj)
e16epi.markers.top5 <- e16epi.markers.top5[order(e16epi.markers.top5$cluster), ]

pdf(file = "e16 epithelium balloon plot.pdf", useDingbats = F, width = 10.5, height = 3.8)
DotPlot(object = e16epithelium, features = c("Epcam",unique(e16epi.markers.top5$gene)), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(angle = 90, size = 14, hjust = 0.5), axis.title = element_blank())
dev.off()


### Generate tSNE plots for Figure
plot1 <- FeaturePlot(e16epithelium, features = c("Krt5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(e16epithelium, features = c("Aqp5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16) )
plot3 <- FeaturePlot(e16epithelium, features = c("Krt14"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-0.4,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(e16epithelium, features = c("Acta2"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes()+ theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-0.4,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("e16 - Krt UMAPS Fig2.pdf", useDingbats = F, width = 5, height = 5)
plot_grid(plotlist = list(plot1,plot2,plot3,plot4), nrow = 2)
dev.off()

plot1 <- FeaturePlot(e16epithelium, features = c("Ung"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(e16epithelium, features = c("Aurka"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16) )
plot3 <- FeaturePlot(e16epithelium, features = c("Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-0.4,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("e16 - Proliferation UMAPS Fig2.pdf", useDingbats = F, width = 5, height = 5)
plot_grid(plotlist = list(plot1,plot2,plot3), nrow = 2)
dev.off()

plot1 <- FeaturePlot(e16epithelium, features = c("Smgc"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16)) 
plot2 <- FeaturePlot(e16epithelium, features = c("Bpifa2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot3 <- FeaturePlot(e16epithelium, features = c("Elf5"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(e16epithelium, features = c("Anxa1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot5 <- FeaturePlot(e16epithelium, features = c("Smoc2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot6 <- FeaturePlot(e16epithelium, features = c("Igfbp2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot7 <- FeaturePlot(e16epithelium, features = c("Tubb5"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot8 <- FeaturePlot(e16epithelium, features = c("Barx1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("e16 - Representative UMAPS Fig2.pdf", useDingbats = F, width = 2.5, height = 17.5)
plot_grid(plotlist = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7,plot8), ncol = 1)
dev.off()

plot1 <- FeaturePlot(e16epithelium, features = c("Igfbp2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,0,0),units = "pt"))
plot2 <- FeaturePlot(e16epithelium, features = c("Smgc"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,0,0),units = "pt"))
plot3 <- FeaturePlot(e16epithelium, features = c("Anxa1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,0,0),units = "pt"))
plot4 <- FeaturePlot(e16epithelium, features = c("Smoc2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,0,0),units = "pt"))
plot5 <- FeaturePlot(e16epithelium, features = c("Bpifa2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,0,0),units = "pt"))
plot6 <- FeaturePlot(e16epithelium, features = c("Elf5"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,0,0),units = "pt"))
plot7 <- FeaturePlot(e16epithelium, features = c("Barx1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,0,0),units = "pt"))
plot8 <- FeaturePlot(e16epithelium, features = c("Tubb5"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,0,0),units = "pt"))
CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7,plot8), ncol = 1)

## visualization of mesenchymal markers
pdf(file = "E12-E16 epithelium (mesenchymal markers).pdf", useDingbats = F, width = 2.5, height = 2.5)
FeaturePlot(e12epithelium, features = "Col1a1", min.cutoff = "q10", max.cutoff = "q90")+ NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-14,0,-0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(e12epithelium, features = "Col3a1", min.cutoff = "q10", max.cutoff = "q90")+ NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-14,0,-0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(e14epithelium, features = "Col1a1", min.cutoff = "q10", max.cutoff = "q90")+ NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-14,0,-0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(e14epithelium, features = "Col3a1", min.cutoff = "q10", max.cutoff = "q90")+ NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-14,0,-0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(e16epithelium, features = "Col1a1", min.cutoff = "q10", max.cutoff = "q90")+ NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-14,0,-0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(e16epithelium, features = "Col3a1", min.cutoff = "q10", max.cutoff = "q90")+ NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-14,0,-0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
dev.off()

## Annotate epithelial clusters
e16epithelium <- SetIdent(e16epithelium, value = "unsupervised")
e16epithelium <- RenameIdents(e16epithelium, '0' = "Myoepithelial", '1' = "End bud",'2'= "Krt19+ duct",'3'= "Basal duct",'4'= "End bud",'5'= "End bud",'6'= "Mesenchyme", '7' = "Myoepithelial", '8'= "Basal duct")
Idents(object = e16epithelium) <- factor(x = Idents(object = e16epithelium), levels = c("End bud", "Krt19+ duct", "Basal duct", "Myoepithelial", "Mesenchyme"))
e16epithelium[["celltype.fixed"]] <- Idents(object = e16epithelium)
e16epithelium <- SetIdent(e16epithelium, value = "celltype.fixed")

e16epithelium <- SetIdent(e16epithelium, value = "unsupervised")
e16epithelium <- RenameIdents(e16epithelium, '0' = "Myoepithelial", '1' = "Smgc+ End bud", '2'= "Krt19+ duct",'3'= "Basal duct",'4'= "Bpifa2+ End bud",'5'= "Proliferative End bud",'6'= "Mesenchyme", '7' = "Proliferative Myoepithelial", '8'= "Proliferative Basal duct")
e16epithelium[["subpopulations"]] <- Idents(object = e16epithelium)
e16epithelium <- SetIdent(e16epithelium, value = "subpopulations")
DimPlot(e16epithelium)

pdf(file = "e16 epithelium re-clustered (Fixed CellType IDs).pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(object = e16epithelium, pt.size = 1.5,label = T,label.size = 5, repel = T, cols = colors, group.by = "celltype.fixed") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank()) + NoLegend()
dev.off()
pdf(file = "e16 epithelium re-clustered (Subpopulation IDs).pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(object = e16epithelium, pt.size = 1.5,label = T,label.size = 5, repel = T, group.by = "subpopulations") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank()) + NoLegend()
dev.off()

pdf(file = "e16 epithelium re-clustered (mesenchymal markers).pdf", useDingbats = F, width = 5.5, height = 5)
FeaturePlot(e16epithelium, features = "Col1a1", min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(e16epithelium, features = "Col3a1", min.cutoff = "q10", max.cutoff = "q90")
dev.off()

e16epithelium <- SetIdent(e16epithelium, value = "celltype.fixed")
e16cell.markers <- FindAllMarkers(e16epithelium, logfc.threshold = 0.25, only.pos = T)
e16cell.markers <- e16cell.markers[e16cell.markers$p_val_adj<0.05, ]
write.csv(e16cell.markers, file = "e16 Cell Type markers (Epithelium).csv")

saveRDS(e16epithelium, file = "E16 epithelium qcd_annotated.rds")




