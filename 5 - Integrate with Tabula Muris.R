# Script for analysis in Figure 9
# ID cells in SMG are integrated with datasets from the tabula muris study

library(dplyr)
library(Seurat)
library(Hmisc)
library(RColorBrewer)
library(viridis)
library(ggplot2)
########### LOAD PREVIOUSLY SAVED FILE WITH EPITHELIAL CLUSTERS FROM ALL STAGES ###########################

load("~/Old 10X projects/GlobalMECAnalysis/Data/droplet_Mammary_seurat_tiss.Robj")
mammary.epithelium <- UpdateSeuratObject(tiss)
load("~/Old 10X projects/GlobalMECAnalysis/Data/droplet_Bladder_seurat_tiss.Robj")
bladder.epithelium <- UpdateSeuratObject(tiss)
load("~/Old 10X projects/GlobalMECAnalysis/Data/droplet_Kidney_seurat_tiss.Robj")
kidney.epithelium <- UpdateSeuratObject(tiss)
load("~/Old 10X projects/GlobalMECAnalysis/Data/droplet_Lung_seurat_tiss.Robj")
lung <- UpdateSeuratObject(tissue.10X)
load("~/Old 10X projects/GlobalMECAnalysis/Data/droplet_Trachea_seurat_tiss.Robj")
trachea <- UpdateSeuratObject(tiss)
load("~/Old 10X projects/GlobalMECAnalysis/Data/facs_Pancreas_seurat_tiss.Robj")
pancreas <- UpdateSeuratObject(tiss)
#load("~/Old 10X projects/GlobalMECAnalysis/Data/droplet_Tongue_seurat_tiss.Robj")
#tongue <- UpdateSeuratObject(tiss)
#tongue[["annotation"]] <- tongue@meta.data$cell_ontology_class
#load("~/Old 10X projects/GlobalMECAnalysis/Data/droplet_Liver_seurat_tiss.Robj")
#liver<- UpdateSeuratObject(tiss)
#liver[["annotation"]] <- liver@meta.data$cell_ontology_class

ID.cells <- readRDS(file = "../../10X revisions/P1-Adult/ID cells subset.rds")
ID.cells[["annotation"]] <- ID.cells$CellType
DimPlot(ID.cells, group.by = "annotation")
ID.cells <- SetIdent(ID.cells, value="annotation")

## rename identities to reflect our findings
ID.cells <- RenameIdents(ID.cells, 'Intercalated duct'='Gfra3 ID', "Smgc+"="Gstt1 ID")
ID.cells[["annotation"]] <- Idents(ID.cells)
## add metadata for SMG tissue
ID.cells[["tissue"]] <- "SMG"

remove(tiss,tissue.10X)

##### Integrate Files to identify similar cell types if present ####
objlist <- list(bladder.epithelium,ID.cells,kidney.epithelium, lung, mammary.epithelium,pancreas,trachea)
objlist <- lapply(X = objlist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

tissue.anchors <- FindIntegrationAnchors(object.list = objlist, dims = 1:30)
tissue.combined <- IntegrateData(anchorset = tissue.anchors, dims = 1:30)
DefaultAssay(tissue.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
tissue.combined <- ScaleData(tissue.combined, verbose = FALSE)
tissue.combined <- RunPCA(tissue.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
tissue.combined <- RunUMAP(tissue.combined, reduction = "pca", dims = 1:30)
tissue.combined <- FindNeighbors(tissue.combined, reduction = "pca", dims = 1:30)

for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,0.9, 1, 1.2, 1.4))
{tissue.combined <- FindClusters(tissue.combined, resolution = resolution)
}
library(clustree)
# Perform non-linear dimensional reduction
pdf('clustree.pdf', width = 30, height = 20)
clustree(tissue.combined, prefix = "integrated_snn_res.")
dev.off() 

tissue.combined <- FindClusters(tissue.combined, resolution = 0.5)
# Visualization
tissue.combined <- SetIdent(tissue.combined, value = "integrated_snn_res.0.5")
pdf("Integration with Tabula muris (res05).pdf", useDingbats = F, width = 6,height = 4)
DimPlot(tissue.combined, group.by = "integrated_snn_res.0.5", label = T, label.size = 5)
dev.off()


tissue.combined <- SetIdent(tissue.combined, value = "annotation")
pdf("Gfra3 and Gstt1 cells (Integration with Tabula muris).pdf", useDingbats = F, width = 6,height = 4)
DimPlot(tissue.combined, group.by = "tissue", label = F) + scale_color_brewer(palette = "Paired")
DimPlot(tissue.combined, reduction = "umap", group.by = "tissue", cells.highlight = WhichCells(tissue.combined, idents = "Gfra3 ID"), cols.highlight = "dodgerblue2")
DimPlot(tissue.combined, reduction = "umap", group.by = "tissue", cells.highlight = WhichCells(tissue.combined, idents = "Gstt1 ID"), cols.highlight = "firebrick2")
dev.off()

tissue.combined <- SetIdent(tissue.combined, value = "tissue")
pdf("Integrated with tabula muris annotated split by tissue.pdf", useDingbats = F, width = 12,height = 12)
DimPlot(subset(tissue.combined, idents = "SMG", invert=T), group.by = "annotation", split.by = "tissue", label = T, repel = T, label.size = 3, ncol = 3) +NoLegend()
dev.off()


## we create new metadata to identify which cluster numbers contain our ID cells
tissue.combined$idcluster <- paste0(tissue.combined$annotation, "_", tissue.combined$integrated_snn_res.0.5)
numberofcells <- as.data.frame(table(tissue.combined$idcluster))
numberofcells <- numberofcells[grep(numberofcells$Var1, pattern = "Gfra3|Gstt1"), ]

library(tidyr)
numberofcells <- numberofcells %>% separate(Var1, c("CellType", "ID", "Cluster"))
library(reshape2)
numberofidcells <- dcast(data = numberofcells,formula = Cluster~CellType,fun.aggregate = sum,value.var = "Freq")

#The resulting table shows that ID cells are in clusters 2 and 6
# We subset those clusters and see what other cells they have
tissue.combined <- SetIdent(tissue.combined, value = "integrated_snn_res.0.5")
clusters.subset <- subset(tissue.combined, idents = c(2,6))

cluster.colors <- c(
  "#c34d54",
  #d7cd9e" 
  #92d8d4"
  "#0076af"
  #6b0081"
)
names(cluster.colors) <- c("2", "6")

pdf("Subset of clusters of interest.pdf", useDingbats = F, width = 4.5, height = 4)
DimPlot(clusters.subset, cols = cluster.colors, pt.size = 1) + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
DimPlot(clusters.subset, pt.size = 1, group.by = "tissue") + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + scale_color_brewer(palette = "Dark2")
dev.off()

## Determine cluster-defining markers for clusters 2 and 6 and visualize top genes
DefaultAssay(tissue.combined) <- "RNA"
unsup.markers <- FindAllMarkers(tissue.combined, logfc.threshold = 0.25, only.pos = T, max.cells.per.ident = 2000, random.seed = 15)

cellmarkerssubset <- unsup.markers[unsup.markers$cluster %in% c(2,6),]
cellmarkerssubset <- cellmarkerssubset[cellmarkerssubset$p_val_adj<0.05, ]
clustermarkers.top <- cellmarkerssubset %>% group_by(cluster) %>% top_n(25,avg_logFC)


library(superheat)
IDgene.matrix <- as.matrix(GetAssayData(subset(tissue.combined, idents = c(2,6))),slot = "data") 
IDgene.matrix <- subset(IDgene.matrix, subset = rownames(IDgene.matrix) %in% unique(clustermarkers.top$gene))

pdf("heatmap clusters 2-6 integrated with TM.pdf", useDingbats = F, width = 5, height = 3)
superheat(t(IDgene.matrix),
          scale = T,legend = T, pretty.order.cols = T, 
          left.label.text.angle = 90,
          bottom.label.text.angle = 90, bottom.label.size = 0.5,  
          bottom.label.text.size =2, left.label.text.size = 4, 
          left.label.col = "white", bottom.label.col = "white",
          grid.hline.size = 0.2, grid.vline.size = 0.2, 
          left.label.size = 0.1, bottom.label.text.alignment = "right", 
          membership.rows = clusters.subset$integrated_snn_res.0.5) 
dev.off()

## heatmap shows clearly that the clusters are heterogeneous and contain subpopulations of cells, likely from the different tissues.
## We next pull the annotations from tabula muris to see which cell types were previously identified in those clusters

clusters.subset[["cluster.annotation"]] <- paste0(clusters.subset$integrated_snn_res.0.5, "_", clusters.subset$annotation, "-", clusters.subset$tissue)
allcellspercluster <- as.data.frame(table(clusters.subset$cluster.annotation))  

allcellspercluster <- allcellspercluster %>% separate(Var1, c("Cluster", "CellType"), sep = "_")
allcellspercluster <- dcast(data = allcellspercluster,formula = Cluster~CellType,fun.aggregate = sum,value.var = "Freq")
allcellspercluster <- allcellspercluster[,-1]
allcellspercluster <- as.data.frame(t(allcellspercluster), stringsAsFactors = T)
names(allcellspercluster) <- c("Cluster 2", "Cluster 6")

write.csv(allcellspercluster, file = "Cells in clusters 2 and 6.csv")


### Annotate plots with the main cell populations in clusters 2 and 6 based on original annotations in Tabula Muris. 
### Only the most abundant populations will be considered

tempdf <- allcellspercluster
tempdf[["celltype"]] <- rownames(tempdf)
tempdf <- tempdf %>% separate(celltype, c("CellType", "Tissue"), sep = "-")

cluster2 <- tempdf[tempdf$`Cluster 2`>0, c(1,3,4)]
cluster2 <- cluster2[order(cluster2$`Cluster 2`, decreasing = T), ][1:10, ]

cluster6 <- tempdf[tempdf$`Cluster 6`>0, c(2:4)]
cluster6 <- cluster6[order(cluster6$`Cluster 6`, decreasing = T), ][1:10, ]

clusters.subset <- SetIdent(clusters.subset, value = "annotation")
pdf("umap split by cluster with cell type annotations (tabula muris).pdf", width = 8.5, height = 4)
DimPlot(subset(clusters.subset, idents = c(cluster2$CellType, cluster6$CellType)), split.by = "integrated_snn_res.0.5") + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
dev.off()

pdf("Cluster 2 and 6 individual gene plots.pdf", useDingbats = F, width = 3, height = 3)
FeaturePlot(tissue.combined, features = c("Gstt1"), min.cutoff = "q10", max.cutoff = "q90") +NoAxes() +NoLegend() + theme(plot.title = element_text(face = "italic", hjust = 0.05))
FeaturePlot(tissue.combined, features = c("Gfra3"), min.cutoff = "q10", max.cutoff = "q90") +NoAxes() +NoLegend()+ theme(plot.title = element_text(face = "italic", hjust = 0.05))
FeaturePlot(tissue.combined, features = c("Cldn4"), min.cutoff = "q10", max.cutoff = "q90") +NoAxes() +NoLegend()+ theme(plot.title = element_text(face = "italic", hjust = 0.05))
FeaturePlot(tissue.combined, features = c("Krt18"), min.cutoff = "q10", max.cutoff = "q90") +NoAxes() +NoLegend()+ theme(plot.title = element_text(face = "italic", hjust = 0.05))
FeaturePlot(tissue.combined, features = c("Krt8"), min.cutoff = "q10", max.cutoff = "q90") +NoAxes() +NoLegend()+ theme(plot.title = element_text(face = "italic", hjust = 0.05))
FeaturePlot(tissue.combined, features = c("Cldn3"), min.cutoff = "q10", max.cutoff = "q90") +NoAxes() +NoLegend()+ theme(plot.title = element_text(face = "italic", hjust = 0.05))
dev.off()




