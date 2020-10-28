# Script for analysis of postnatal SMG epithelium 
# associated with Figure 3-4
# This scrip also contains code for nalysis of ID cells in Figure 8

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
            "#fbe54b",
            "#e74134",
            "hotpink")

cell.types <- c("Acinar", "Ascl3+ duct", "Basal duct", "Bpifa2+", "Bpifa2+ Proacinar", "End bud", "Endothelial",
                "Erythroid", "GCT", "Glial cells", "Intercalated duct", "Krt19+ duct", "Macrophages", "Mast cells",
                "Mesenchyme", "Mitotic cells", "Myoepithelial", "Nerves", "NK cells", "Smgc+", "Smgc+ Proacinar",
                "Smooth muscle", "Striated duct", "Stromal", "Smgc+ Male", "Smgc+ Female")

names(colors) <- cell.types
colors

stage.colors <- c(
  "#c34d54",
  "#d7cd9e",
  "#019f84",
  "#92d8d4",
  "#0076af",
  "#6b0081"
)
names(stage.colors) <- c("p1", "p30", "adult", "P1", "P30", "Adult")

########################## Load previously annotated SEURAT files############################
#############################################################################################

p.epithelium <- readRDS("../Postnatal epithelium.rds")
DefaultAssay(p.epithelium) <- "RNA"

pdf(file = "postnatal epithelium (subset).pdf", useDingbats = F, width = 6, height = 4)
DimPlot(object = p.epithelium, pt.size = 0.5,label =F,repel = F, label.size = 5, group.by = "CellType", cols = colors) + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank())
dev.off()

### Split object to re-cluster epithelium from individual developmental stages
p.split <- SplitObject(p.epithelium, split.by = "stage")
### standard worflow for normalization and scaling
for (i in 1:length(p.split)) {
 p.split[[i]] <- NormalizeData(object =p.split[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
 p.split[[i]] <- FindVariableFeatures(object =p.split[[i]], selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x =p.split[[i]])
 p.split[[i]] <- ScaleData(object =p.split[[i]], features = all.genes, vars.to.regress = "percent.mt")
 p.split[[i]] <- RunPCA(object =p.split[[i]], features = VariableFeatures(object =p.split[[i]]))
}

# extract individual stages after normalization and scaling
p1epithelium <-p.split[[1]]
p30epithelium <-p.split[[2]]
adultepithelium <-p.split[[3]]

save(p1epithelium, p30epithelium, adultepithelium, file = "Postnatal epithelium re-scaled before reclustering (backup).RData")

######## ANALYSIS OF p1 epithelium
#Evaluate optimal clustering resolution using clustree package
DefaultAssay(p1epithelium) <- "RNA"
ElbowPlot(object = p1epithelium)
p1epithelium <- FindNeighbors(object = p1epithelium, dims = 1:15,force.recalc = T)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,0.9, 1, 1.2, 1.4)){
  p1epithelium <- FindClusters(p1epithelium, resolution = resolution)
}
pdf('p1_clustree.pdf', width = 15, height = 10)
clustree(p1epithelium,prefix = "RNA_snn_res.")
dev.off()
# optimal resolution at 0.6.
# perform clustering
p1epithelium <- FindClusters(object = p1epithelium, resolution = 0.6)
p1epithelium <- RunUMAP(object = p1epithelium, dims = 1:15)

pdf(file = "p1 epithelium re-clustered (original label).pdf", useDingbats = F, width = 5.5, height = 4)
DimPlot(object = p1epithelium, pt.size = 1,label = T,repel = F, label.size = 5, group.by = "CellType", cols = colors) + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank())
dev.off()
pdf(file = "p1 epithelium re-clustered (new cluster IDs).pdf", useDingbats = F, width = 3, height = 3)
DimPlot(object = p1epithelium, pt.size = 1,label = T,repel = F, label.size = 5, group.by = "RNA_snn_res.0.5") + NoAxes() + NoLegend()
dev.off()


p1epithelium$unsupervised <- Idents(p1epithelium)

### Identify p1 markers for putative epithelial subpopulations
p1epithelium <- SetIdent(p1epithelium, value = "unsupervised")
p1epi.markers <- FindAllMarkers(p1epithelium, logfc.threshold = 0.25, only.pos = T)
write.csv(p1epi.markers, file = "p1 Unsupervised Epithelium markers.csv")

#Re-order identities for easier visualization of subclusters from similar populations
target <- c(1,7,8,0,3,6,5,2,4,9)
Idents(p1epithelium) <- factor(Idents(p1epithelium), levels = target)
DimPlot(p1epithelium)

p1epi.markers$cluster <- factor(p1epi.markers$cluster, levels = target)
p1epi.markers.top5 <- p1epi.markers %>% group_by(cluster) %>% top_n(-5,p_val_adj)
p1epi.markers.top5 <- p1epi.markers.top5[order(p1epi.markers.top5$cluster),]

pdf(file = "p1 epithelium balloon plot.pdf", useDingbats = F, width = 11, height = 3)
DotPlot(object = p1epithelium, features = c("Epcam",unique(p1epi.markers.top5$gene)), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"), axis.text.y = element_text(angle = 90, size = 14, hjust = 0.5), axis.title = element_blank()) + NoLegend()
dev.off()

### Generate UMAP plots for Figure
plot1 <- FeaturePlot(p1epithelium, features = c("Bpifa2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16)) 
plot2 <- FeaturePlot(p1epithelium, features = c("Dcpp1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot3 <- FeaturePlot(p1epithelium, features = c("Smgc"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(p1epithelium, features = c("Ssr4"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot5 <- FeaturePlot(p1epithelium, features = c("Wfdc18"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot6 <- FeaturePlot(p1epithelium, features = c("Ccna2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot7 <- FeaturePlot(p1epithelium, features = c("Krt7"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot8 <- FeaturePlot(p1epithelium, features = c("Fbln2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot9 <- FeaturePlot(p1epithelium, features = c("Cnn1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("P1 Representative UMAPS Fig.pdf", useDingbats = F, width = 2.5, height = 20)
plot_grid(plotlist = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7,plot8,plot9), ncol = 1)
dev.off()

plot1 <- FeaturePlot(p1epithelium, features = c("Ung"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-10,0,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(p1epithelium, features = c("Aurka"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-10,0,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16) )
plot3 <- FeaturePlot(p1epithelium, features = c("Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-10,0,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("P1 - Proliferation UMAPS Fig.pdf", useDingbats = F, width = 5, height = 5)
plot_grid(plotlist = list(plot1,plot2,plot3), nrow = 2)
dev.off()

plot1 <- FeaturePlot(p1epithelium, features = c("Aqp5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(p1epithelium, features = c("Bhlha15"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot3 <- FeaturePlot(p1epithelium, features = c("Gstt1"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(p1epithelium, features = c("Dcdc2a"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("P1 - F7 UMAPS-1.pdf", useDingbats = F, width = 2, height = 8)
plot_grid(plotlist = list(plot1,plot2,plot3,plot4), ncol =1)
dev.off()

plot1 <- FeaturePlot(p1epithelium, features = c("Gfra3"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(p1epithelium, features = c("Kit"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot3 <- FeaturePlot(p1epithelium, features = c("Nkd2"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(p1epithelium, features = c("Esp18"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("P1 - F7 UMAPS-2.pdf", useDingbats = F, width = 2, height = 8)
plot_grid(plotlist = list(plot1,plot2,plot3,plot4), ncol =1)
dev.off()

### Annotate epithelial clusters
p1epithelium <- SetIdent(p1epithelium, value = "unsupervised")
p1epithelium <- RenameIdents(p1epithelium, '0' = "Smgc+ Proacinar", '1' = "Bpifa2+ Proacinar",'2'= "Krt19+ duct",'3'= "Smgc+ Proacinar",
                             '4'= "Basal duct",'5'= "Mitotic cells",'6'= "Smgc+ Proacinar", '7' = "Bpifa2+ Proacinar", '8'= "Bpifa2+ Proacinar",
                             '9'= "Myoepithelial")
Idents(object = p1epithelium) <- factor(x = Idents(object = p1epithelium), levels = c("Bpifa2+ Proacinar", "Smgc+ Proacinar", "Krt19+ duct", "Basal duct", "Myoepithelial", "Mitotic cells"))
p1epithelium[["celltype.fixed"]] <- Idents(object = p1epithelium)
p1epithelium <- SetIdent(p1epithelium, value = "celltype.fixed")
p1epithelium[["subpopulations"]] <- Idents(object = p1epithelium)
p1epithelium <- SetIdent(p1epithelium, value = "subpopulations")
DimPlot(p1epithelium)

pdf(file = "p1 epithelium re-clustered (Fixed CellType IDs).pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(object = p1epithelium, pt.size = 1.5,label = T,label.size = 5, repel = T, cols = colors, group.by = "celltype.fixed") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank()) + NoLegend()
dev.off()
pdf(file = "p1 epithelium re-clustered (Subpopulation IDs).pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(object = p1epithelium, pt.size = 1.5,label = T,label.size = 5, repel = T, group.by = "subpopulations") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank()) + NoLegend()
dev.off()

pdf(file = "p1 epithelium re-clustered (mesenchymal markers).pdf", useDingbats = F, width = 5.5, height = 5)
FeaturePlot(p1epithelium, features = "Col1a1", min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(p1epithelium, features = "Col3a1", min.cutoff = "q10", max.cutoff = "q90")
dev.off()

p1epithelium <- SetIdent(p1epithelium, value = "celltype.fixed")
p1cell.markers <- FindAllMarkers(p1epithelium, logfc.threshold = 0.25, only.pos = T)
p1cell.markers <-p1cell.markers[p1cell.markers$p_val_adj <0.05, ]
write.csv(p1cell.markers, file = "P1 Epithelial cell type markers.csv")


saveRDS(p1epithelium, file = "p1 epithelium qcd_annotated.rds")

### ANALYSIS OF p30 epithelium
### Evaluate optimal clustering resolution using clustree package
DefaultAssay(p30epithelium) <- "RNA"
ElbowPlot(object = p30epithelium)
p30epithelium <- FindNeighbors(object = p30epithelium, dims = 1:12,force.recalc = T)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,0.9, 1, 1.2, 1.4)){
  p30epithelium <- FindClusters(p30epithelium, resolution = resolution)
}
pdf('p30_clustree.pdf', width = 15, height = 10)
clustree(p30epithelium,prefix = "RNA_snn_res.")
dev.off()

# perform clustering WITH RESOLUTION OF 0.9 == lower resolution did not properly discriminate between basal duct and GCT cells.
# myoepithelial cells and basal duct cells clustering together regardless of resolution likely due to low number of cells
p30epithelium <- FindClusters(object = p30epithelium, resolution = 0.9)
p30epithelium <- RunUMAP(object = p30epithelium, dims = 1:12)
p30epithelium$unsupervised <- Idents(p30epithelium)

pdf(file = "p30 epithelium re-clustered (original label).pdf", useDingbats = F, width = 5.5, height = 4)
DimPlot(object = p30epithelium, pt.size = 1,label = T,repel = T, label.size = 5, group.by = "CellType", cols = colors) + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank())
dev.off()

pdf(file = "p30 epithelium re-clustered (colored by sex).pdf", useDingbats = F, width = 5.5, height = 4)
DimPlot(object = p30epithelium, pt.size = 1,group.by = "sex", cols = c("#00009980", "#FD5353")) + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank())
dev.off()


pdf(file = "p30 epithelium re-clustered (new cluster IDs).pdf", useDingbats = F, width = 3, height = 3)
DimPlot(object = p30epithelium, pt.size = 1,label = T,repel = F, label.size = 5, group.by = "unsupervised") +NoAxes()  +NoLegend()
dev.off()


### Identify p30 markers for putative epithelial subpopulations
p30epithelium <- SetIdent(p30epithelium, value = "unsupervised")
p30epi.markers <- FindAllMarkers(p30epithelium, logfc.threshold = 0.25, only.pos = T)
write.csv(p30epi.markers, file = "p30 Unsupervised Epithelium markers.csv")

#Re-order identities for easier visualization of subclusters from similar populations
target <- c(0,2,7,1,4,6,5,3,9,8)
Idents(p30epithelium) <- factor(Idents(p30epithelium), levels = target)
DimPlot(p30epithelium)

p30epi.markers$cluster <- factor(p30epi.markers$cluster, levels =  target)
p30epi.markers.top5 <- p30epi.markers %>% group_by(cluster) %>% top_n(-5,p_val_adj)
p30epi.markers.top5 <- p30epi.markers.top5[order(p30epi.markers.top5$cluster), ]

pdf(file = "p30 epithelium balloon plot2.pdf", useDingbats = F, width = 11, height = 3)
DotPlot(object = p30epithelium, features = c("Epcam",unique(p30epi.markers.top5$gene)), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"), axis.text.y = element_text(angle = 90, size = 14, hjust = 0.5), axis.title = element_blank())
dev.off()

### Generate UMAP plots for Figure
plot1 <- FeaturePlot(p30epithelium, features = c("Lpo"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16)) 
plot2 <- FeaturePlot(p30epithelium, features = c("Mucl2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot3 <- FeaturePlot(p30epithelium, features = c("Smgc"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(p30epithelium, features = c("Serpinb11"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot5 <- FeaturePlot(p30epithelium, features = c("Gfra3"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot6 <- FeaturePlot(p30epithelium, features = c("Ngf"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot7 <- FeaturePlot(p30epithelium, features = c("Ldhb"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot8 <- FeaturePlot(p30epithelium, features = c("Fxyd2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot9 <- FeaturePlot(p30epithelium, features = c("Ascl3"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot10 <- FeaturePlot(p30epithelium, features = c("Krt5"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("P30 Representative UMAPS Fig.pdf", useDingbats = F, width = 2.5, height = 20)
plot_grid(plotlist = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7,plot8,plot9,plot10), ncol = 1)
dev.off()


plot1 <- FeaturePlot(p30epithelium, features = c("Ung"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-10,0,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(p30epithelium, features = c("Aurka"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-10,0,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16) )
plot3 <- FeaturePlot(p30epithelium, features = c("Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-10,0,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("p30 - Proliferation UMAPS Fig.pdf", useDingbats = F, width = 5, height = 5)
plot_grid(plotlist = list(plot1,plot2,plot3), nrow = 2)
dev.off()

plot1 <- FeaturePlot(p30epithelium, features = c("Aqp5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin= unit(c(-10,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(p30epithelium, features = c("Bhlha15"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-10,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot3 <- FeaturePlot(p30epithelium, features = c("Smgc"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-10,0,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(p30epithelium, features = c("Dcdc2a"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-10,0,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("p30 - F7 UMAPS-1.pdf", useDingbats = F, width = 2, height = 8)
plot_grid(plotlist = list(plot1,plot2,plot3,plot4), ncol =1, align = "hv")
dev.off()

plot1 <- FeaturePlot(p30epithelium, features = c("Gfra3"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(p30epithelium, features = c("Kit"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot3 <- FeaturePlot(p30epithelium, features = c("Nkd2"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(p30epithelium, features = c("Esp18"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("p30 - F7 UMAPS-2.pdf", useDingbats = F, width = 2, height = 8)
plot_grid(plotlist = list(plot1,plot2,plot3,plot4), ncol =1)
dev.off()

## Annotate epithelial clusters
p30epithelium <- SetIdent(p30epithelium, value = "unsupervised")
DimPlot(p30epithelium, split.by = "sex") # this allows us to see that clusters 1 and 4 are similar cell types but different in males and females
p30epithelium <- RenameIdents(p30epithelium, '0' = "Acinar", '1' = "Smgc+ Female",'2'= "Acinar",'3'= "Striated duct",'4'= "Smgc+ Male",'5'= "GCT",'6'= "Intercalated duct",
                              "7" = "Acinar", "8" = "Basal duct", '9' = "Ascl3+ duct")
#label MECs manually because they clustered together with basal duct cells
MECs <- WhichCells(p30epithelium, expression = Cnn1>0&Acta2>0)
Idents(p30epithelium, cells= MECs) <- "Myoepithelial"
DimPlot(p30epithelium)
Idents(object = p30epithelium) <- factor(x = Idents(object = p30epithelium), levels = c("Acinar", "Smgc+ Female", "Smgc+ Male", "Intercalated duct", "Striated duct", "GCT", "Ascl3+ duct", "Basal duct", "Myoepithelial"))
p30epithelium[["celltype.fixed"]] <- Idents(object = p30epithelium)
p30epithelium <- SetIdent(p30epithelium, value = "celltype.fixed")

p30epithelium <- SetIdent(p30epithelium, value = "unsupervised")
p30epithelium <- RenameIdents(p30epithelium, '0' = "Acinar 1", '1' = "Smgc+ Female",'2'= "Acinar 2",'3'= "Striated duct",'4'= "Smgc+ Male",'5'= "GCT",'6'= "Intercalated duct",
                              "7" = "Acinar 3", "8" = "Basal duct", '9' = "Ascl3+ duct")
Idents(p30epithelium, cells= MECs) <- "Myoepithelial"
Idents(object = p30epithelium) <- factor(x = Idents(object = p30epithelium), levels = c("Acinar 1","Acinar 2", "Acinar 3", "Smgc+ Female", "Smgc+ Male", "Intercalated duct", "Striated duct", "GCT", "Ascl3+ duct", "Basal duct", "Myoepithelial"))
p30epithelium[["subpopulations"]] <- Idents(object = p30epithelium)
p30epithelium <- SetIdent(p30epithelium, value = "subpopulations")
DimPlot(p30epithelium)

## The above analysis shows a sexually dimorphic population that expresses Smgc in females and Serpinb1 in males
pdf("Smgc expression at P30 - Epithelial clusters only - split by sex.pdf", useDingbats = F, width = 9, height = 4.5) 
FeaturePlot(p30epithelium, features = c("Smgc"), min.cutoff = "q10", max.cutoff = "q90", split.by = "sex") 
FeaturePlot(p30epithelium, features = c("Serpinb11"), min.cutoff = "q10", max.cutoff = "q90", split.by = "sex") 
FeaturePlot(p30epithelium, features = c("Gstt1"), min.cutoff = "q10", max.cutoff = "q90", split.by = "sex")
dev.off()

pdf("GCT markers at P30 - Epithelial clusters only - split by sex.pdf", useDingbats = F, width = 9, height = 4.5) 
FeaturePlot(p30epithelium, features = c("Ngf"), min.cutoff = "q10", max.cutoff = "q90", split.by = "sex") 
FeaturePlot(p30epithelium, features = c("Egf"), min.cutoff = "q10", max.cutoff = "q90", split.by = "sex") 
dev.off()

pdf("Sexual dimorphic markers at P30 - Epithelial clusters only - split by sex.pdf", useDingbats = F, width = 9, height = 4.5) 
FeaturePlot(p30epithelium, features = c("Dcdc2a"), min.cutoff = "q10", max.cutoff = "q90", split.by = "sex") 
FeaturePlot(p30epithelium, features = c("Cdkn1c"), min.cutoff = "q10", max.cutoff = "q90", split.by = "sex") 
FeaturePlot(p30epithelium, features = c("Gstt1"), min.cutoff = "q10", max.cutoff = "q90", split.by = "sex") 
FeaturePlot(p30epithelium, features = c("Serpinb11"), min.cutoff = "q10", max.cutoff = "q90", split.by = "sex") 
dev.off()

pdf(file = "p30 epithelium re-clustered (Fixed CellType IDs).pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(object = p30epithelium, pt.size = 1.5,label = T,label.size = 5, repel = T, cols = colors, group.by = "celltype.fixed") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank()) + NoLegend()
dev.off()
pdf(file = "p30 epithelium re-clustered (Subpopulation IDs).pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(object = p30epithelium, pt.size = 1.5,label = T,label.size = 5, repel = T, group.by = "subpopulations") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank()) + NoLegend()
dev.off()


## make a metadata slot with fixed cell identities per sex for differential expression analysis between males and females
p30epithelium <- SetIdent(p30epithelium, value = "celltype.fixed")
p30epithelium[["temp.label"]] <- p30epithelium$celltype.fixed
p30epithelium <- SetIdent(p30epithelium, value = "temp.label")
p30epithelium <-RenameIdents(p30epithelium, "Smgc+ Female" = "Smgc+", "Smgc+ Male" = "Smgc+")
p30epithelium$temp.label <- Idents(p30epithelium)
DimPlot(p30epithelium)

p30epithelium[["celltype.fixed.sex"]] <- paste0(p30epithelium$temp.label, "_", p30epithelium$sex)
p30epithelium <- SetIdent(p30epithelium, value = "celltype.fixed.sex")
write.csv(as.data.frame(table(Idents(p30epithelium))), file = "P30 epithelium _ fixed cell counts.csv")

p30epithelium <- SetIdent(p30epithelium, value = "celltype.fixed")
p30cell.markers <- FindAllMarkers(p30epithelium, logfc.threshold = 0.25, only.pos = T)
p30cell.markers <-p30cell.markers[p30cell.markers$p_val_adj <0.05, ]
p30cell.markers.top <- p30cell.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
write.csv(p30cell.markers, file = "p30 Epithelial cell type markers.csv")

saveRDS(p30epithelium, file = "p30 epithelium qcd_annotated.rds")

### IDENTIFY DIFFERENTIALLY EXPRESSED GENES BETWEEN MALES AND FEMALES ###

DefaultAssay(p30epithelium) <- "RNA"
p30epithelium <- SetIdent(p30epithelium, value = "celltype.fixed.sex")
#order identities alphabetically to facilitate differential expression analysis with a for loop
Idents(p30epithelium) <- factor(Idents(p30epithelium), levels = sort(levels(Idents(p30epithelium)),decreasing = F))
DimPlot(p30epithelium)
p30epithelium$celltype.fixed.sex <- Idents(p30epithelium)

### Differential expression analysis
comparisons.list <- list()
j=0
for (i in seq(from = 1,to =  (length(levels(Idents(p30epithelium)))), by = 2)) { 
  comparisons.list[[j+1]] <- FindMarkers(p30epithelium, ident.1 = levels(Idents(p30epithelium))[i], ident.2 = levels(Idents(p30epithelium))[i+1], only.pos = F, logfc.threshold = 0.5, min.pct = 0.25)  #only genes present in at least 25% of cells
  comparisons.list[[j+1]]$gene <- rownames(comparisons.list[[j+1]])
  names(comparisons.list[[j+1]])[3] <- levels(Idents(p30epithelium))[i]
  names(comparisons.list[[j+1]])[4] <- levels(Idents(p30epithelium))[i+1]
  comparisons.list[[j+1]]<-comparisons.list[[j+1]][comparisons.list[[j+1]]$p_val_adj <0.05, ]
  comparisons.list[[j+1]]<-comparisons.list[[j+1]][order(comparisons.list[[j+1]]$avg_logFC), ]
  j=j+1
}
## create new sorting order (alphabetically to ensure coherence between exported files)
p30epithelium <- SetIdent(p30epithelium, value = "temp.label")
Idents(p30epithelium) <- factor(Idents(p30epithelium), levels = sort(levels(Idents(p30epithelium)),decreasing = F))

#### Export tables with sex-dependent genes
for (i in 1:length(comparisons.list)) {
  identity <- levels(Idents(p30epithelium))[[i]]
  write.csv(x = comparisons.list[[i]],file = paste0("../P1-Adult/MvsF gene lists/",identity, " M vs F.csv"),row.names = T)
}
#### Determine number of differentially expressed genes per cell type
number.of.degs.per.cell <- c()
numDEGs.mvsf <- plyr::ldply(comparisons.list, rbind)
numDEGs.mvsf <- numDEGs.mvsf[,-c(1:2, 5:6)]
numDEGs.mvsf <- numDEGs.mvsf[,-c(1,3,5,7,9,11,13,15,17,19,21,23)]
numDEGs.mvsf <- as.data.frame(sapply(numDEGs.mvsf, function(x) sum(!(is.na(x)))))
names(numDEGs.mvsf)[1] <- "degs"
rownames(numDEGs.mvsf)<-sub("_Male", "", rownames(numDEGs.mvsf))
names(comparisons.list) <- rownames(numDEGs.mvsf) # add list names for later
write.csv(numDEGs.mvsf, "Number of DEGs per cell type between males and females.csv")

pdf("Bar graph of sex-dependent DEGs in P30 SMG.pdf", useDingbats = F, width = 5,height = 5)
par(mar=c(12,4,2,2))
barplot(numDEGs.mvsf$degs,names.arg = rownames(numDEGs.mvsf),
        ylab="Number of genes", main = "Sex-dependent genes per cell type",
        border="black", col = "blue", las=2, axis.lty = 1)
dev.off()


## The following code exports violin plots for the top up and down-regulated genes
## per cell type (if less than 10 are present, only those are shown)
plot.grid <- list()
p30epithelium <- SetIdent(p30epithelium, value = "temp.label")
Idents(p30epithelium) <- factor(Idents(p30epithelium), levels = sort(levels(Idents(p30epithelium)),decreasing = F))
DimPlot(p30epithelium)
##for cycle with number of identities (8)
for (i in 1:8) {
  identity <- levels(Idents(p30epithelium))[[i]]
  #male genes
  comparisons.list[[i]]<-comparisons.list[[i]][order(comparisons.list[[i]]$avg_logFC), ]
  plot.list <- list()
  if(sum(comparisons.list[[i]]$avg_logFC<0)>0){
    for (j in 1:sum(comparisons.list[[i]]$avg_logFC<0)) {
      plot.list[[j]] <- VlnPlot(p30epithelium, idents = identity, features = comparisons.list[[i]]$gene[j],pt.size = 0,cols = c("#000099", "#FF6666"),split.by = "sex") +NoLegend() +theme(axis.title = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(size = 10,face = "italic")) +NoAxes()
    }
    if(length(plot.list)>10){
      ggsave(filename = paste0(identity," male genes.pdf"),dpi = 300, units = "in", width = 20,height = 1.5,plot = plot_grid(plotlist = plot.list[1:10],ncol = 20)) 
    }
    if(length(plot.list)<=10){
      ggsave2(filename = paste0(identity," male genes.pdf"),dpi = 300, units = "in", width = 20,height = 1.5,plot = plot_grid(plotlist = plot.list,ncol = 20)) 
    }
  }  
  #female genes
  comparisons.list[[i]]<-comparisons.list[[i]][order(comparisons.list[[i]]$avg_logFC,decreasing = T), ]
  plot.list <- list()
  if(sum(comparisons.list[[i]]$avg_logFC>0)>0){
    for (j in 1:sum(comparisons.list[[i]]$avg_logFC>0)) {
      plot.list[[j]] <- VlnPlot(p30epithelium, idents = identity, features = comparisons.list[[i]]$gene[j],pt.size = 0,cols = c("#000099", "#FF6666"),split.by = "sex") +NoLegend() +theme(axis.title = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(size = 10,face = "italic")) +NoAxes()
    }
    if(length(plot.list)>10){
      ggsave2(filename = paste0(identity," female genes.pdf"),dpi = 300, units = "in", width = 20,height = 1.5,plot = plot_grid(plotlist = plot.list[1:10],ncol = 20)) 
    }
    if(length(plot.list)<=10){
      ggsave2(filename = paste0(identity," female genes.pdf"),dpi = 300, units = "in", width = 20,height = 1.5,plot = plot_grid(plotlist = plot.list,ncol = 20)) 
    }
  }}


######## Analysis of adult epithelium
#Evaluate optimal clustering resolution using clustree package
DefaultAssay(adultepithelium) <- "RNA"
ElbowPlot(object = adultepithelium)
adultepithelium <- FindNeighbors(object = adultepithelium, dims = 1:12,force.recalc = T)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,0.9, 1, 1.2, 1.4)){
  adultepithelium <- FindClusters(adultepithelium, resolution = resolution)
}
pdf('adult_clustree.pdf', width = 15, height = 10)
clustree(adultepithelium,prefix = "RNA_snn_res.")
dev.off()
# optimal resolution is 0.5 although striated duct and GCTs are not discriminated
# If resolution is increased, artifacts occur in other populations.


# perform clustering
adultepithelium <- FindClusters(object = adultepithelium, resolution = 0.5)
adultepithelium <- RunUMAP(object = adultepithelium, dims = 1:12)
adultepithelium$unsupervised <- Idents(adultepithelium)

pdf(file = "adult epithelium re-clustered (original label).pdf", useDingbats = F, width = 5.5, height = 4)
DimPlot(object = adultepithelium, pt.size = 1,label = T,repel = T, label.size = 5, group.by = "CellType", cols = colors) + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank())
dev.off()
pdf(file = "adult epithelium re-clustered (new cluster IDs).pdf", useDingbats = F, width = 3, height = 3)
DimPlot(object = adultepithelium, pt.size = 1,label = T,repel = F, label.size = 5, group.by = "unsupervised") + NoAxes() + NoLegend()
dev.off()

### Identify adult markers for putative epithelial subpopulations
adultepithelium <- SetIdent(adultepithelium, value = "unsupervised")
adultepi.markers <- FindAllMarkers(adultepithelium, logfc.threshold = 0.25, only.pos = T)
write.csv(adultepi.markers, file = "adult Unsupervised Epithelium markers.csv")

#Re-order identities for easier visualization of subclusters from similar populations
target <- c(0,4,7,3,10,1,5,2,9,6,11,8)
Idents(adultepithelium) <- factor(Idents(adultepithelium), levels = target)
DimPlot(adultepithelium)

adultepi.markers$cluster <- factor(adultepi.markers$cluster, levels =  target)
adultepi.markers.top5 <- adultepi.markers %>% group_by(cluster) %>% top_n(-5,p_val_adj)
adultepi.markers.top5 <- adultepi.markers.top5[order(adultepi.markers.top5$cluster), ]

pdf(file = "adult epithelium balloon plot.pdf", useDingbats = F, width = 12, height =4)
DotPlot(object = adultepithelium, features = c("Epcam",unique(adultepi.markers.top5$gene)), cols = "Spectral", col.min = 0, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90,face = "italic", hjust = 1, vjust = 0.5), axis.text.y = element_text(angle = 90, size = 14, hjust = 0.5), axis.title = element_blank()) + NoLegend()
dev.off()


### Generate tSNE plots for Figure
plot1 <- FeaturePlot(adultepithelium, features = c("Krt5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(adultepithelium, features = c("Aqp5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16) )
plot3 <- FeaturePlot(adultepithelium, features = c("Krt14"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-0.4,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(adultepithelium, features = c("Acta2"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes()+ theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-0.4,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("adult - Krt UMAPS Fig.pdf", useDingbats = F, width = 5, height = 5)
plot_grid(plotlist = list(plot1,plot2,plot3,plot4), nrow = 2)
dev.off()

plot1 <- FeaturePlot(adultepithelium, features = c("Ung"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(adultepithelium, features = c("Aurka"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,0,-0.4,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16) )
plot3 <- FeaturePlot(adultepithelium, features = c("Mki67"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-0.4,0,0,0),units = "cm"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("adult - Proliferation UMAPS Fig.pdf", useDingbats = F, width = 5, height = 5)
plot_grid(plotlist = list(plot1,plot2,plot3), nrow = 2)
dev.off()

plot1 <- FeaturePlot(adultepithelium, features = c("Prol1"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16)) 
plot2 <- FeaturePlot(adultepithelium, features = c("Bpifa2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot3 <- FeaturePlot(adultepithelium, features = c("Muc19"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(adultepithelium, features = c("Smgc"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot5 <- FeaturePlot(adultepithelium, features = c("Gfra3"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot6 <- FeaturePlot(adultepithelium, features = c("Fxyd2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot7 <- FeaturePlot(adultepithelium, features = c("Ascl3"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot8 <- FeaturePlot(adultepithelium, features = c("Smoc2"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot9 <- FeaturePlot(adultepithelium, features = c("Ifi202b"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot10 <- FeaturePlot(adultepithelium, features = c("Myl9"), pt.size = 0.5,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-18,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("adult - Representative UMAPS Fig.pdf", useDingbats = F, width = 2.5, height = 20)
plot_grid(plotlist = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7,plot8, plot9, plot10), ncol = 1)
dev.off()

plot1 <- FeaturePlot(adultepithelium, features = c("Aqp5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(adultepithelium, features = c("Bhlha15"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot3 <- FeaturePlot(adultepithelium, features = c("Gstt1"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(adultepithelium, features = c("Dcdc2a"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("adult - F7 UMAPS-1.pdf", useDingbats = F, width = 2, height = 8)
plot_grid(plotlist = list(plot1,plot2,plot3,plot4), ncol =1)
dev.off()

plot1 <- FeaturePlot(adultepithelium, features = c("Gfra3"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot2 <- FeaturePlot(adultepithelium, features = c("Kit"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot3 <- FeaturePlot(adultepithelium, features = c("Nkd2"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
plot4 <- FeaturePlot(adultepithelium, features = c("Esp18"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-16,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
pdf("adult - F7 UMAPS-2.pdf", useDingbats = F, width = 2, height = 8)
plot_grid(plotlist = list(plot1,plot2,plot3,plot4), ncol =1)
dev.off()

## Annotate epithelial clusters
# Transfer GCT annotation because it was not discriminated from SD cells
adultepithelium <- SetIdent(adultepithelium, value = "CellType")
GCTcells <- WhichCells(adultepithelium, idents = "GCT")
                       
adultepithelium <- SetIdent(adultepithelium, value = "unsupervised")
adultepithelium <- RenameIdents(adultepithelium, '0' = "Acinar", '1' = "Smgc+",'2'= "Striated duct",'3'= "Bpifa2+",
                                '4'= "Acinar",'5'= "Intercalated duct",'6'= "Basal duct", '7' = "Acinar", '8'= "Myoepithelial",
                                '9'= "Ascl3+ duct",'10'= "Bpifa2+",'11'= "Basal duct")
Idents(adultepithelium, cells =GCTcells) <- "GCT"

Idents(object = adultepithelium) <- factor(x = Idents(object = adultepithelium), levels = c("Acinar", "Bpifa2+", "Smgc+", "Intercalated duct", "Striated duct", "GCT", "Ascl3+ duct", "Basal duct", "Myoepithelial"))
adultepithelium[["celltype.fixed"]] <- Idents(object = adultepithelium)

adultepithelium <- SetIdent(adultepithelium, value = "unsupervised")
adultepithelium <- RenameIdents(adultepithelium, '0' = "Acinar", '1' = "Smgc+",'2'= "Striated duct",'3'= "Bpifa2+",
                                '4'= "Acinar 2",'5'= "Intercalated duct",'6'= "Basal duct 1", '7' = "Acinar 3", '8'= "Myoepithelial",
                                '9'= "Ascl3+ duct",'10'= "Muc19+",'11'= "Basal duct 2")
Idents(adultepithelium, cells =GCTcells) <- "GCT"
Idents(object = adultepithelium) <- factor(x = Idents(object = adultepithelium), levels = c("Acinar", "Acinar 2", "Acinar 3", "Bpifa2+", "Muc19+", "Smgc+", "Intercalated duct", "Striated duct","GCT", "Ascl3+ duct", "Basal duct 1", "Basal duct 2", "Myoepithelial"))
adultepithelium[["subpopulations"]] <- Idents(object = adultepithelium)
adultepithelium <- SetIdent(adultepithelium, value = "subpopulations")
DimPlot(adultepithelium)

pdf(file = "adult epithelium re-clustered (Fixed CellType IDs).pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(object = adultepithelium, pt.size = 1.5,label = T,label.size = 5, repel = T, cols = colors, group.by = "celltype.fixed") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank()) + NoLegend()
dev.off()
pdf(file = "adult epithelium re-clustered (Subpopulation IDs).pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(object = adultepithelium, pt.size = 1.5,label = T,label.size = 5, repel = T, group.by = "subpopulations") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_blank()) + NoLegend()
dev.off()

adultepithelium <- SetIdent(adultepithelium, value = "celltype.fixed")
adultcell.markers <- FindAllMarkers(adultepithelium, logfc.threshold = 0.25, only.pos = T)
adultcell.markers <-adultcell.markers[adultcell.markers$p_val_adj <0.05, ]
adultcell.markers.top <- adultcell.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
write.csv(adultcell.markers, file = "adult Epithelial cell type markers.csv")

saveRDS(p30epithelium, file = "p30 epithelium qcd_annotated.rds")

saveRDS(adultepithelium, file = "adult epithelium qcd_annotated.rds")

pdf("Dot plots with top cell markers in P30 and adult.pdf", useDingbats = F, width = 10, height = 4)
DotPlot(p30epithelium, features = unique(p30cell.markers.top$gene), dot.scale = 5,group.by = "celltype.fixed") + theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())
DotPlot(adultepithelium, features = unique(adultcell.markers.top$gene), dot.scale = 5,group.by = "celltype.fixed") + theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5), axis.title = element_blank())
dev.off()


## visualization of mesenchymal markers
pdf(file = "P1-adult epithelium (mesenchymal markers).pdf", useDingbats = F, width = 2.5, height = 2.5)
FeaturePlot(p1epithelium, features = "Col1a1", min.cutoff = "q10", max.cutoff = "q90")+ NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-14,0,-0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(p1epithelium, features = "Col3a1", min.cutoff = "q10", max.cutoff = "q90")+ NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-14,0,-0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(p30epithelium, features = "Col1a1", min.cutoff = "q10", max.cutoff = "q90")+ NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-14,0,-0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(p30epithelium, features = "Col3a1", min.cutoff = "q10", max.cutoff = "q90")+ NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-14,0,-0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(adultepithelium, features = "Col1a1", min.cutoff = "q10", max.cutoff = "q90")+ NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-14,0,-0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(adultepithelium, features = "Col3a1", min.cutoff = "q10", max.cutoff = "q90")+ NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(-14,0,-0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
dev.off()

####################################################
### Comparison of proacinar and Acinar cells ###

p.epithelium <- readRDS("../Postnatal epithelium.rds")
DefaultAssay(p.epithelium) <- "RNA"
DimPlot(p.epithelium)

acinarmarkers <- p30cell.markers[p30cell.markers$cluster %in% "Acinar", ]
acinarmarkers <- acinarmarkers[order(acinarmarkers$avg_logFC, decreasing = T), ]

pdf("Proacinar - Acinar shared markers.pdf", useDingbats = F, width = 8, height = 2.2)
DotPlot(p.epithelium, features = acinarmarkers$gene, 
        idents = c("Acinar", "Bpifa2+ Proacinar", "Smgc+ Proacinar"), dot.min = 0.05, dot.scale = 4) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"), axis.title = element_blank())
dev.off( )


#### Markers of acinar subtypes (serous vs mucous)

pdf("UMAP with final annotations - adult gland.pdf", useDingbats = F, width = 6,height = 5)
DimPlot(adultepithelium, group.by = "celltype.fixed", cols = colors)
dev.off()

pdf("Serous-Mucous markers.pdf",useDingbats = F, width = 3,height = 3)
FeaturePlot(adultepithelium, features = "Prol1", min.cutoff = "q10", max.cutoff = "q90") + NoAxes() + theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(adultepithelium, features = "Mucl2", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(adultepithelium, features = "Dcpp1", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(adultepithelium, features = "Dcpp2", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(adultepithelium, features = "Bpifa2", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(adultepithelium, features = "Car6", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(adultepithelium, features = "Lpo", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(adultepithelium, features = "Cd44", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(adultepithelium, features = "Alcam", min.cutoff = "q10", max.cutoff = "q90")+ NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
dev.off()

pdf("Serous-Mucous markers_p30.pdf",useDingbats = F, width = 3,height = 3)
FeaturePlot(p30epithelium, features = "Prol1", min.cutoff = "q10", max.cutoff = "q90") + NoAxes() + theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(p30epithelium, features = "Mucl2", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(p30epithelium, features = "Dcpp1", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(p30epithelium, features = "Dcpp2", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(p30epithelium, features = "Bpifa2", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(p30epithelium, features = "Car6", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(p30epithelium, features = "Lpo", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(p30epithelium, features = "Cd44", min.cutoff = "q10", max.cutoff = "q90") + NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
FeaturePlot(p30epithelium, features = "Alcam", min.cutoff = "q10", max.cutoff = "q90")+ NoAxes()+ theme(plot.title = element_text(face="italic", hjust = 0.05))
dev.off()

####################################################
### Analysis of GFRA3 and SMGC Intercalated duct ###
## Figure 8

# identify differentially expressed markers between Smgc+ ID and GFRa3+ ID
gfra3.vs.smgc <- FindMarkers(p.epithelium, ident.1 = "Smgc+", ident.2 = "Intercalated duct", only.pos = F,logfc.threshold = 0.25,min.pct = 0.25)
gfra3.vs.smgc <- gfra3.vs.smgc[gfra3.vs.smgc$p_val_adj<0.05, ]
gfra3.vs.smgc <- gfra3.vs.smgc[order(gfra3.vs.smgc$avg_logFC), ]
write.csv(gfra3.vs.smgc, file = "Gfra3 vs Smgc cells.csv")

#select top differentially expressed genes and create expression matrix for visualization in heatmap
genes.to.plot <- rbind(gfra3.vs.smgc %>% top_n(-15, avg_logFC), gfra3.vs.smgc %>% top_n(15, avg_logFC) )
AvgExp.all.cells <- AverageExpression(p.epithelium)
AvgExpRNA <- AvgExp.all.cells$RNA
AvgExpRNA <- AvgExpRNA[rownames(AvgExpRNA) %in% rownames(genes.to.plot), colnames(AvgExpRNA) %in% c("Smgc+", "Intercalated duct", "Smgc+ Proacinar", "Krt19+ duct")]

pdf(file = "Smgc Gfra3 heatmap.pdf", useDingbats = F, width = 6.5, height = 4.5)
superheat(t(AvgExpRNA),
         scale = T,legend = T, pretty.order.cols = T,
         left.label.text.angle = 90,
         bottom.label.text.angle = 90, bottom.label.size = 0.5,  
         bottom.label.text.size = 4, left.label.text.size = 4, 
         left.label.col = "white", bottom.label.col = "white",
         grid.hline.size = 0.2, grid.vline.size = 0.2, 
         left.label.size = 0.1, bottom.label.text.alignment = "right") 
dev.off()


gfra3.features <- FindMarkers(p.epithelium, ident.1 = "Intercalated duct", group.by = "CellType", only.pos = T, min.pct = 0.25)
gfra3.features$gene <- rownames(gfra3.features)
gfra3.features<-gfra3.features[gfra3.features$p_val_adj<0.05, ]
gfra3.features.top100 <- gfra3.features[1:100,]
write.csv(gfra3.features, file = "Gfra3-ID genes in postnatal epithelium.csv")

smgc.features <- FindMarkers(p.epithelium, ident.1 = "Smgc+", group.by = "CellType", only.pos = T,min.pct = 0.25)
smgc.features$gene <- rownames(smgc.features)
smgc.features<-smgc.features[smgc.features$p_val_adj<0.05, ]
smgc.features.top100 <- smgc.features[1:100,]
write.csv(smgc.features, file = "Gstt1-ID genes in postnatal epithelium.csv")

allIDgenes <- unique(c(gfra3.features$gene, smgc.features$gene)) # get a list of all defining genes for ID populations and remove duplicates
table(Idents(p.epithelium)) #there are 804 cells in smgc/gstt1 and ID clusters
IDgene.matrix <- as.matrix(GetAssayData(subset(p.epithelium, idents = c("Intercalated duct", "Smgc+")),slot = "counts")) # resulting matrix should have 804 columns
IDgene.matrix <- subset(IDgene.matrix, subset = rownames(IDgene.matrix) %in% allIDgenes)

superheat(t(IDgene.matrix),
          scale = T,legend = T, pretty.order.cols = T, 
          left.label.text.angle = 90,
          bottom.label.text.angle = 90, bottom.label.size = 0.5,  
          bottom.label.text.size = 4, left.label.text.size = 4, 
          left.label.col = "white", bottom.label.col = "white",
          grid.hline.size = 0.2, grid.vline.size = 0.2, 
          left.label.size = 0.1, bottom.label.text.alignment = "right") 


#### SEPARATE SMGC AND ID CELLS TO MERGE WITH TABULA MURIS FILES ####
## For figures 8-9 ##
DimPlot(p.epithelium)
ID.cells <- subset(p.epithelium, idents = c("Smgc+", "Intercalated duct"))
saveRDS(ID.cells, file = "ID cells subset.rds")


pdf("Figure plots.pdf", useDingbats = F, width = 2, height = 2)
FeaturePlot(p1epithelium, features = c("Aqp5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(p1epithelium, features = c("Bhlha15"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(p1epithelium, features = c("Smgc"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(p1epithelium, features = c("Dcdc2a"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(p30epithelium, features = c("Aqp5"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90") + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(p30epithelium, features = c("Bhlha15"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(p30epithelium, features = c("Smgc"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
FeaturePlot(p30epithelium, features = c("Dcdc2a"), pt.size = 1,min.cutoff = "q10", max.cutoff = "q90")  + NoLegend() +NoAxes() + theme(title = element_blank(), panel.border = element_rect(size = 0.5,colour = "black"), plot.margin=unit(c(0,3,0,0),units = "pt"), plot.title = element_text(face = "italic", hjust = 0.02, vjust = -8, size = 16))
dev.off()


