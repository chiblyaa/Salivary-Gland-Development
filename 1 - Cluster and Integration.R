# This script will import expression data into SEURAT and will perform 
# analysis for Figure 1 in manuscript

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(Hmisc)
library(cowplot)

## Load gene expression data
## Expression matrix files available at GSE150327
e12smg.data <- Read10X(data.dir = "~/E12_V3/outs/filtered_gene_bc_matrices/mm10plus/")
e14smg.data <- Read10X(data.dir = "~/E14_V3/outs/filtered_gene_bc_matrices/mm10plus/")
e16smg.data <- Read10X(data.dir = "~/SalGland_E16_v2/filtered_gene_bc_matrices/mm10plus/")
p1smg.data <- Read10X(data.dir = "~/SalGland_P1_v3/filtered_gene_bc_matrices/mm10plus/")
p30.male.data <- Read10X_h5(filename =  "~/P30_M-selected/P30_M_filtered_gene_bc_matrices_h5.h5")
p30.female.data <- Read10X_h5(filename = "~/P30_F-selected/P30_F_filtered_gene_bc_matrices_h5.h5")
adult.data <-  Read10X(data.dir = "~/Adult SMG/filtered_feature_bc_matrix/")

# Define a color scheme for each cell type to generate consistent plots 
# We will use 24 distinct colors to annotate different cell populations based on identified cell types:
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

cell.type.levels <- c("End bud", "Krt19+ duct", "Basal duct", "Myoepithelial", "Bpifa2+ Proacinar", "Smgc+ Proacinar", "Mitotic cells",
                      "Bpifa2+", "Acinar", "Smgc+", "Intercalated duct", "Ascl3+ duct", "GCT", "Striated duct",    
                      "Endothelial", "Smooth muscle",
                      "Nerves",  "Glial cells",   "Macrophages", "Mast cells", "NK cells","Mesenchyme", "Stromal", "Erythroid") # use this for consistent sorting of cell types in plots (optional)

#############################################################################################
########################## ANALYSIS OF E12 EPITHELIUM #######################################
#############################################################################################

e12smg <- CreateSeuratObject(counts = e12smg.data, min.cells = 3, min.features = 200)
e12smg[["percent.mt"]] <- PercentageFeatureSet(object = e12smg, pattern = "^mt-")
VlnPlot(object = e12smg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = e12smg, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = e12smg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
e12smg <- subset(x = e12smg, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 5)
e12smg <- NormalizeData(object = e12smg, normalization.method = "LogNormalize", scale.factor = 10000)
e12smg <- FindVariableFeatures(object = e12smg, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = e12smg), 10)
plot1 <- VariableFeaturePlot(object = e12smg)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(x = e12smg)
e12smg <- ScaleData(object = e12smg, features = all.genes)
e12smg <- RunPCA(object = e12smg, features = VariableFeatures(object = e12smg))
ElbowPlot(object = e12smg)
e12smg <- FindNeighbors(object = e12smg, dims = 1:12)
e12smg <- FindClusters(object = e12smg, resolution = 0.6)
e12smg <- RunTSNE(object = e12smg, dims = 1:12)
e12smg <- RunUMAP(e12smg, dims = 1:12)

pdf("E12 SMG UMAP unsupervised.pdf", width = 3.5, height = 3, useDingbats = F)
DimPlot(object = e12smg, reduction = "umap", pt.size = 1,label = T,label.size = 5, repel = F, group.by = "seurat_clusters") + NoLegend()
dev.off()

#### Identify cluster markers for annotation
e12smg <- SetIdent(e12smg, value = "seurat_clusters")
e12markers <- FindAllMarkers(e12smg, only.pos = T,logfc.threshold = 0.25)
e12markers.top5 <- e12markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.csv(e12markers, file = "E12 SMG Unsupervised markers (SEURAT).csv")

DotPlot(e12smg, features = unique(e12markers.top5$gene), cols = "Spectral", dot.scale = 4,group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 9), axis.title = element_blank(), axis.text.y = element_text(size = 9)) + NoLegend()

pdf("Balloon plot E12 annotation strategy.pdf", useDingbats = F, width = 4,height = 3.2)
DotPlot(e12smg, features = c("Epcam", "Krt14", "Krt5", "Sox9", "Krt19", "Aqp5", "Acta2", "Cnn1", "Pecam1", "Tubb3", "Ncam1", "Vim", "Col1a1", "Twist1", "Adgre1", "Alas2"), dot.scale = 4, dot.min = 0.05, group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 9), axis.title = element_blank(), axis.text.y = element_text(size = 9)) + NoLegend()
dev.off()

e12smg[["stage"]] <- "E12"
e12smg <- RenameIdents(e12smg, '0' = "Mesenchyme",'1' = "Mesenchyme",'2' = "Mesenchyme",'3' = "Mesenchyme",'4' = "Mesenchyme",'5' = "Mesenchyme",'6' = "Mesenchyme", 
                       '7'= "End bud", '8' = "Krt19+ duct", '9'= "Nerves", '10'="Smooth muscle", '11' = "Erythroid", '12'="Endothelial", '13'="Macrophages")
Idents(e12smg) <- factor(Idents(e12smg), levels = sort(levels(Idents(e12smg)),decreasing = F))
DimPlot(e12smg)
e12smg[["CellType"]] <- Idents(e12smg)

pdf("E12 SMG UMAP annotated.pdf", width = 3.5, height = 3, useDingbats = F)
DimPlot(e12smg, group.by = "CellType", label = F, label.size = 5, pt.size = 0.5, repel = T, cols = colors) + theme(axis.text = element_text(size=14), axis.title = element_blank()) + NoLegend()
dev.off()

e12smg <- SetIdent(e12smg, value = "CellType")
e12cellcounts <- as.data.frame(table(Idents(e12smg)))
write.csv(e12cellcounts, file = "E12 cell counts.csv")
e12.cell.markers <- FindAllMarkers(e12smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
e12.cell.markers.top5 <- e12.cell.markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.csv(e12.cell.markers, file = "E12 SMG Cell Type markers (SEURAT).csv")
DotPlot(e12smg, features = unique(e12.cell.markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "CellType") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

saveRDS(e12smg, file = "E12 SMG Annotated (SEURAT v3).rds")

############################################################################################
########################## ANALYSIS OF E14 EPITHELIUM #######################################
#############################################################################################
e14smg <- CreateSeuratObject(counts = e14smg.data, min.cells = 3, min.features = 200)
e14smg[["percent.mt"]] <- PercentageFeatureSet(object = e14smg, pattern = "^mt-")
VlnPlot(object = e14smg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = e14smg, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = e14smg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
e14smg <- subset(x = e14smg, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 5)
e14smg <- NormalizeData(object = e14smg, normalization.method = "LogNormalize", scale.factor = 10000)
e14smg <- FindVariableFeatures(object = e14smg, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = e14smg), 10)
plot1 <- VariableFeaturePlot(object = e14smg)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(x = e14smg)
e14smg <- ScaleData(object = e14smg, features = all.genes)
e14smg <- RunPCA(object = e14smg, features = VariableFeatures(object = e14smg))
ElbowPlot(object = e14smg)
e14smg <- FindNeighbors(object = e14smg, dims = 1:12)
e14smg <- FindClusters(object = e14smg, resolution = 0.9)
e14smg <- RunTSNE(object = e14smg, dims = 1:12)
e14smg <- RunUMAP(e14smg, dims = 1:12)
DimPlot(object = e14smg, reduction = "tsne", pt.size = 1,label = T,label.size = 6)

pdf("E14 SMG UMAP unsupervised.pdf", width = 3.5, height = 3, useDingbats = F)
DimPlot(object = e14smg, reduction = "umap", pt.size = 1,label = T,label.size = 5, repel = T, group.by = "seurat_clusters") +NoLegend()
dev.off()

e14smg <- SetIdent(e14smg, value = "seurat_clusters") # make sure identity is set to unsupervised clusters
e14markers <- FindAllMarkers(e14smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
e14markers.top5 <- e14markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.csv(e14markers, file = "e14 SMG Unsupervised markers (SEURAT).csv")
DotPlot(e14smg, features = unique(e14markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

pdf("Balloon plot E14 annotation strategy.pdf", useDingbats = F, width = 4,height = 3.2)
DotPlot(e14smg, features = c("Epcam", "Krt14", "Krt5", "Sox9", "Krt19", "Acta2", "Cnn1", "Pecam1", "Tubb3", "Ncam1", "Vim", "Col1a1", "Twist1", "Adgre1", "Alas2"),dot.min = 0.05, dot.scale = 4,group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 9), axis.title = element_blank(), axis.text.y = element_text(size = 9)) + NoLegend()
dev.off()


e14smg[["stage"]] <- "E14"
e14smg <- RenameIdents(e14smg, '0' = "Mesenchyme",'1' = "Mesenchyme",'3' = "Mesenchyme",'14' = "Mesenchyme",'8' = "Mesenchyme",
                       '19' = "Endothelial",'16' = "Nerves",'17'="Macrophages",'18'="Macrophages",'11' = "Erythroid",'15' = "Erythroid", '13'="Glial cells",
                       '7'= "Krt19+ duct", '6' = "End bud", '5'="End bud", '2'="End bud",'9'="End bud",'4'="End bud",'12'="End bud",'10'="Basal duct")
Idents(e14smg) <- factor(Idents(e14smg), levels = sort(levels(Idents(e14smg)),decreasing = F))
DimPlot(e14smg)
e14smg[["CellType"]] <- Idents(e14smg)

pdf("E14 SMG UMAP annotated.pdf", width = 5, height = 3.2, useDingbats = F)
DimPlot(e14smg, group.by = "CellType", label = F, label.size = 5, pt.size = 0.5, repel = T, cols = colors) + theme(axis.text = element_text(size=14), axis.title = element_blank())
dev.off()

e14smg <- SetIdent(e14smg, value = "CellType")
e14cellcounts <- as.data.frame(table(Idents(e14smg)))
write.csv(e14cellcounts, file = "E14 cell counts.csv")
e14.cell.markers <- FindAllMarkers(e14smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
e14.cell.markers.top5 <- e14.cell.markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.csv(e14.cell.markers, file = "e14 SMG Cell Type markers (SEURAT).csv")
DotPlot(e14smg, features = unique(e14.cell.markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "CellType") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

saveRDS(e14smg, file = "E14 SMG Annotated (SEURAT v3).rds")

#############################################################################################
########################## ANALYSIS OF E16 EPITHELIUM #######################################
#############################################################################################

e16smg <- CreateSeuratObject(counts = e16smg.data, min.cells = 3, min.features = 200)
e16smg[["percent.mt"]] <- PercentageFeatureSet(object = e16smg, pattern = "^mt-")
VlnPlot(object = e16smg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = e16smg, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = e16smg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
e16smg <- subset(x = e16smg, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 5)
e16smg <- NormalizeData(object = e16smg, normalization.method = "LogNormalize", scale.factor = 10000)
e16smg <- FindVariableFeatures(object = e16smg, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = e16smg), 10)
plot1 <- VariableFeaturePlot(object = e16smg)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(x = e16smg)
e16smg <- ScaleData(object = e16smg, features = all.genes)
e16smg <- RunPCA(object = e16smg, features = VariableFeatures(object = e16smg))
ElbowPlot(object = e16smg)
e16smg <- FindNeighbors(object = e16smg, dims = 1:12)
e16smg <- FindClusters(object = e16smg, resolution = 0.9)
e16smg <- RunTSNE(object = e16smg, dims = 1:12)
e16smg <- RunUMAP(object = e16smg, dims = 1:12)
DimPlot(object = e16smg, reduction = "tsne", pt.size = 1,label = T,label.size = 6)

pdf("E16 SMG UMAP unsupervised.pdf", width = 3.5, height = 3, useDingbats = F)
DimPlot(object = e16smg, reduction = "umap", pt.size = 1,label = T,label.size = 5, repel = F, group.by = "seurat_clusters") +NoLegend()
dev.off()

e16smg <- SetIdent(e16smg, value = "seurat_clusters") # make sure identity is set to unsupervised clusters
e16markers <- FindAllMarkers(e16smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
write.csv(e16markers, file = "E16 SMG Unsupervised markers (SEURAT).csv", sep = "\t", row.names = F)
e16markers.top5 <- e16markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
DotPlot(e16smg, features = unique(e16markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

pdf("Balloon plot E16 annotation strategy.pdf", useDingbats = F, width = 4,height = 3.2)
DotPlot(e16smg, features = c("Epcam", "Krt14", "Krt5", "Sox9", "Krt19", "Acta2", "Cnn1", "Pecam1", "Tubb3", "Ncam1", "Vim", "Col1a1", "Twist1", "Adgre1", "Alas2"), dot.min = 0.05,   dot.scale = 4,group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 9), axis.title = element_blank(), axis.text.y = element_text(size = 9)) + NoLegend()
dev.off()

e16smg[["stage"]] <- "E16"
e16smg <- RenameIdents(e16smg, '0' = "Mesenchyme",'1' = "Mesenchyme",'2' = "Mesenchyme",'3' = "Mesenchyme",'6' = "Mesenchyme",'8' = "Mesenchyme",'15' = "Mesenchyme",'16' = "Mesenchyme",
                       '12' = "Erythroid", '11'="Macrophages",'18'="Macrophages", '21' = "Nerves", '20'="Glial cells", '19' = "Mast cells", '17'="Smooth muscle", '14'="Basal duct", '13'="Krt19+ duct",
                       '10'= "Endothelial", '9' = "Myoepithelial", '5'="Erythroid", '4'="End bud",'7'="Erythroid")
Idents(e16smg) <- factor(Idents(e16smg), levels = sort(levels(Idents(e16smg)),decreasing = F))
DimPlot(e16smg)
e16smg[["CellType"]] <- Idents(e16smg)

pdf("E16 SMG UMAP annotated.pdf", width = 5, height = 3.2, useDingbats = F)
DimPlot(e16smg, group.by = "CellType", label = F, label.size = 5, pt.size = 0.5, cols = colors, repel = T) + theme(axis.text = element_text(size=14), axis.title = element_blank())
dev.off()

e16smg <- SetIdent(e16smg, value = "CellType")
e16cellcounts <- as.data.frame(table(Idents(e16smg)))
write.csv(e16cellcounts, file = "E16 cell counts.csv")
e16.cell.markers <- FindAllMarkers(e16smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
e16.cell.markers.top5 <- e16.cell.markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.csv(e16.cell.markers, file = "e16 SMG Cell Type markers (SEURAT).csv")
DotPlot(e16smg, features = unique(e16.cell.markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "CellType") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

saveRDS(e16smg, file = "E16 SMG Annotated (SEURAT v3).rds")

#############################################################################################
########################## ANALYSIS OF p1 EPITHELIUM #######################################
#############################################################################################

p1smg <- CreateSeuratObject(counts = p1smg.data, min.cells = 3, min.features = 200)
p1smg[["percent.mt"]] <- PercentageFeatureSet(object = p1smg, pattern = "^mt-")
VlnPlot(object = p1smg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = p1smg, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = p1smg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
p1smg <- subset(x = p1smg, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 5)
p1smg <- NormalizeData(object = p1smg, normalization.method = "LogNormalize", scale.factor = 10000)
p1smg <- FindVariableFeatures(object = p1smg, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = p1smg), 10)
plot1 <- VariableFeaturePlot(object = p1smg)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(x = p1smg)
p1smg <- ScaleData(object = p1smg, features = all.genes)
p1smg <- RunPCA(object = p1smg, features = VariableFeatures(object = p1smg))
ElbowPlot(object = p1smg)
p1smg <- FindNeighbors(object = p1smg, dims = 1:12)
p1smg <- FindClusters(object = p1smg, resolution = 0.9)
p1smg <- RunTSNE(object = p1smg, dims = 1:12)
DimPlot(object = p1smg, reduction = "tsne", pt.size = 1,label = T,label.size = 6)
p1smg <- RunUMAP(object = p1smg, dims = 1:12)

pdf("P1 SMG UMAP unsupervised.pdf", width = 3.5, height = 3, useDingbats = F)
DimPlot(object = p1smg, reduction = "umap", pt.size = 0.5,label = T,label.size = 5, repel = F, group.by = "seurat_clusters") + NoLegend()
dev.off()

p1smg <- SetIdent(p1smg, value = "seurat_clusters") # make sure identity is set to unsupervised clusters
p1markers <- FindAllMarkers(p1smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
write.csv(p1markers, file = "P1 SMG Unsupervised markers (SEURAT).csv")
p1markers.top5 <- p1markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
DotPlot(p1smg, features = unique(p1markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

pdf("Balloon plot P1 annotation strategy.pdf", useDingbats = F, width = 5,height = 3.2)
DotPlot(p1smg, features = c("Epcam", "Smgc", "Aqp5", "Bhlha15", "Krt19", "Klk1", "Cftr", "Ascl3", "Krt14", "Krt5", "Acta2", "Cnn1", "Pecam1", "Tubb3", "Ncam1", "Vim", "Col1a1", "Twist1", "Adgre1", "Alas2"),dot.min = 0.05, dot.scale = 4,group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 9), axis.title = element_blank(), axis.text.y = element_text(size = 9)) + NoLegend()
dev.off()

p1smg[["stage"]] <- "P1"
p1smg <- RenameIdents(p1smg, '0' = "Bpifa2+ Proacinar",'9' = "Bpifa2+ Proacinar",'10' = "Bpifa2+ Proacinar",'1' = "Smgc+ Proacinar", '2' = "Smgc+ Proacinar",'3' = "Smgc+ Proacinar",
                      '17'= "Endothelial", '5'="Krt19+ duct", '6' = "Smgc+ Proacinar", '7' = "Basal duct", '8'= "Mitotic cells", '11' = "Bpifa2+ Proacinar", '12'="Krt19+ duct",
                      '15' = "Basal duct", '14' = "Myoepithelial", '16'="Macrophages", '4' = "Mesenchyme", '13'="Mesenchyme", '18' ="Erythroid")
Idents(p1smg) <- factor(Idents(p1smg), levels = sort(levels(Idents(p1smg)),decreasing = F))
DimPlot(p1smg)
p1smg[["CellType"]] <- Idents(p1smg)

pdf("P1 SMG UMAP Annotated.pdf", width = 5, height = 3.2, useDingbats = F)
DimPlot(p1smg, group.by = "CellType", label = F, label.size = 5, pt.size = 1, repel = T, cols = colors) + theme(axis.text = element_text(size=14), axis.title = element_blank())
dev.off()

p1smg <- SetIdent(p1smg, value = "CellType")
p1cellcounts <- as.data.frame(table(Idents(p1smg)))
write.csv(p1cellcounts, file = "P1 cell counts.csv")
p1.cell.markers <- FindAllMarkers(p1smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
p1.cell.markers.top5 <- p1.cell.markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.csv(p1.cell.markers, file = "P1 SMG Cell Type markers (SEURAT).csv")
DotPlot(p1smg, features = unique(p1.cell.markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "CellType") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

saveRDS(p1smg, file = "P1 SMG Annotated (SEURAT v3).rds")


#################################################################################################
########################## Analysis of P30 and adults ###########################################
#################################################################################################

## Goal: Because similar cell types are expected between P30 and adults, these datasets are integrated 
## to evaluate batch effects and to designate consistent cell type annotations in all three datasets. 

p30male <- CreateSeuratObject(p30.male.data,min.cells = 3, min.features = 200)
p30female <- CreateSeuratObject(p30.female.data, min.cells = 3, min.features = 200)
adultsmg <- CreateSeuratObject(counts = adult.data, min.cells = 3, min.features = 200)
p30male[["sex"]] <- "Male"
p30female[["sex"]] <- "Female"
adultsmg[["sex"]] <-"Female"
p30male[["stage"]] <- "P30"
p30female[["stage"]] <- "P30"
adultsmg[["stage"]] <- "Adult"
p30male[["sample"]] <- "S1"
p30female[["sample"]] <- "S2"
adultsmg[["sample"]] <- "S3"

smgP.integrated.list <- list(p30male, p30female, adultsmg)

for (i in 1:length(smgP.integrated.list)) {
  smgP.integrated.list[[i]][["percent.mt"]] <- PercentageFeatureSet(smgP.integrated.list[[i]], pattern = "^mt-")
  smgP.integrated.list[[i]] <- subset(smgP.integrated.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 10)
  smgP.integrated.list[[i]] <- NormalizeData(object = smgP.integrated.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  smgP.integrated.list[[i]] <- FindVariableFeatures(smgP.integrated.list[[i]], selection.method = "vst", nfeatures = 2000)
}
smgP.integrated.anchors <- FindIntegrationAnchors(smgP.integrated.list, dims = 1:30)
smgP.integrated <- IntegrateData(smgP.integrated.anchors, dims = 1:30)

remove(smgP.integrated.list)

# Run the standard workflow for visualization and clustering
smgP.integrated <- ScaleData(smgP.integrated, verbose = FALSE)
smgP.integrated <- RunPCA(smgP.integrated, npcs = 30)
smgP.integrated <- RunUMAP(smgP.integrated, reduction = "pca", dims = 1:30)
smgP.integrated <- FindNeighbors(smgP.integrated, reduction = "pca", dims = 1:30)
smgP.integrated <- FindClusters(smgP.integrated, resolution = 0.9)

pdf(file = "UMAP P30_Adult Integrated_unsupervised.pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(smgP.integrated, reduction = "umap", label = T, label.size = 5, repel = F, group.by = "seurat_clusters", split.by = "sex") + NoLegend()
dev.off()

pdf("Balloon plot P30-Adult annotation strategy.pdf", useDingbats = F, width = 5,height = 3.2)
DotPlot(smgP.integrated, features = c("Epcam", "Smgc", "Aqp5", "Bhlha15", "Krt19", "Klk1", "Cftr", "Ascl3", "Krt14", "Krt5", "Acta2", "Cnn1", "Pecam1", "Tubb3", "Ncam1", "Vim", "Col1a1", "Twist1", "Adgre1", "Alas2"), dot.min = 0.05, dot.scale = 4,group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 9), axis.title = element_blank(), axis.text.y = element_text(size = 9)) + NoLegend()
dev.off()


smgP.integrated <- SetIdent(smgP.integrated, value = "seurat_clusters")
smgP.integrated.int.markers <- FindAllMarkers(smgP.integrated, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
smgP.integrated <- RenameIdents(smgP.integrated, '0' = "Smgc+", '1' = "Acinar", '2' = "Acinar",
                                '3' = "Striated duct", '4' = "Macrophages", '5' = "Endothelial",
                                '6'= "Intercalated duct", '7'="Bpifa2+", '8' = "Striated duct",
                                '9' = "Basal duct", '10'= "Myoepithelial", '11'="GCT",
                                '12'="Ascl3+ duct", '13'="NK cells", '14'="Macrophages", '15'="Erythroid",
                                '16'="Bpifa2+", '17'="Basal duct")
write.csv(smgP.integrated.int.markers, file = "P30_Adult SMG unsupervised markers.csv")

DefaultAssay(smgP.integrated) <- "RNA"

#a small population of stromal cells was identified clustered together with myoepithelial cells
# they are manually separated based on expression of known stromal markers 
stromalcells <- WhichCells(smgP.integrated, expression = Vim>0&Col3a1>0, slot = "counts") # select Stromal cells based on expression of Collagen and Vimentin
Idents(smgP.integrated, cells=stromalcells) <- "Stromal" 

DimPlot(smgP.integrated)
Idents(smgP.integrated) <- factor(Idents(smgP.integrated), levels = sort(levels(Idents(smgP.integrated)),decreasing = F))
smgP.integrated[["CellType"]] <- Idents(smgP.integrated)

pdf(file = "UMAP P30_Adult Integrated_annotated.pdf", useDingbats = F, width = 5, height = 3.2)
DimPlot(smgP.integrated, group.by = "CellType", label = F, label.size = 5, pt.size = 0.5, repel = T, cols = colors) + theme(axis.text = element_text(size=14), axis.title = element_blank())
dev.off()

pdf(file = "UMAP P30_Adult split by stage.pdf", useDingbats = F, width = 10, height = 5)
DimPlot(smgP.integrated, reduction = "umap", split.by = "stage", pt.size = 0.5, group.by= "stage", cols = c("navy", "hotpink2")) + NoLegend() + theme(axis.title = element_blank(), axis.text = element_text(size = 12))
dev.off()

smgP.integrated <- SetIdent(smgP.integrated, value = "CellType")
smgP.integrated[["celltype.stage"]] <- paste0(smgP.integrated$stage, "_", smgP.integrated$CellType)
smgP.integrated[["celltype.stage.sex"]] <- paste0(smgP.integrated$stage, "_", smgP.integrated$sex, "_", smgP.integrated$CellType)

smgP.integratedcounts <- as.data.frame(table(smgP.integrated$celltype.stage))
smgP.integratedcounts2 <- as.data.frame(table(smgP.integrated$celltype.stage.sex))
write.csv(smgP.integratedcounts2, file = "Number of cells per cluster in P30 and Adults.csv")

saveRDS(smgP.integrated, file = "P30_adult integrated Annotated (SEURAT v3).rds")


### Separate P30 from adults for visualization and to generate list of defining genes
smgP.integrated <-SetIdent(smgP.integrated, value = "stage")
p30subset <- subset(smgP.integrated, idents = "P30")
p30subset <- SetIdent(p30subset, value = "CellType")

pdf("P30 SMG UMAP unsupervised.pdf", useDingbats = F, height = 5, width = 5.5)
DimPlot(p30subset, group.by = "seurat_clusters", pt.size = 1, label = T, label.size = 5, repel = F) + theme(axis.text = element_text(size = 12)) + NoLegend()
dev.off()
pdf("P30 SMG UMAP colored by sex.pdf",  useDingbats = F, height = 5, width = 5.5)
DimPlot(p30subset, label = F, group.by = "sex", pt.size = 0.5, cols = c("#00009980", "#FD5353")) + theme(axis.title = element_blank(), axis.text = element_text(size = 12))
dev.off()
pdf("P30 SMG UMAP colored by cell Type.pdf",  useDingbats = F, height = 3, width = 3.5)
DimPlot(p30subset, pt.size = 0.5, label = F, label.size = 5, repel = T, cols = colors, group.by = "CellType") + theme(axis.title = element_blank(), axis.text = element_text(size = 12)) + NoLegend()
dev.off()
pdf("Smgc Expression in P30 males and females (UMAP).pdf", useDingbats = F, height = 4.5, width = 10)
FeaturePlot(p30subset, features = "Smgc", split.by = "sex", ncol = 1, min.cutoff = "q10")
dev.off()
p30subset[["CellType.sex"]] <- paste0(p30subset@meta.data$CellType, "_", p30subset@meta.data$sex)
saveRDS(p30subset, "P30_Male_and_female combined - annotated (split from integrated).rds")

p30subset <- SetIdent(p30subset, value="CellType")
p30.cell.markers <- FindAllMarkers(p30subset, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
write.csv(p30.cell.markers, file = "P30 SMG Cell Type markers (SEURAT).csv")

adultsubset <- subset(smgP.integrated, idents = "Adult")
adultsubset <- SetIdent(adultsubset, value="CellType")
adult.cell.markers <- FindAllMarkers(adultsubset, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
write.csv(adult.cell.markers, file = "Adult SMG Cell Type markers (SEURAT).csv")


##### Extract individual male and female SEURAT objects with proper cell annotations
##### UMAPs will also be generated for individual stages to upload to SGMAP
tempseurat <- SplitObject(smgP.integrated, split.by = "sample")
#S1 = p30 male
#s2 = p30 female
#s3 = adult

p30male <- tempseurat$S1
p30female <- tempseurat$S2
adultsmg <- tempseurat$S3
remove(tempseurat)

pdf("P30 male UMAP annotated.pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(p30male, pt.size = 1, label = T, label.size = 5, repel = T, cols = colors, group.by = "CellType") + theme(axis.title = element_blank(), axis.text = element_text(size = 12)) + NoLegend()
dev.off()

pdf("P30 female UMAP annotated.pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(p30female, pt.size = 1, label = T, label.size = 5, repel = T, cols = colors, group.by = "CellType") + theme(axis.title = element_blank(), axis.text = element_text(size = 12)) + NoLegend()
dev.off()

pdf("Adult SMG UMAP unsupervised.pdf", useDingbats = F, width = 5.5, height = 5)
DimPlot(adultsmg, group.by = "seurat_clusters", label = T, label.size = 5, repel = F, pt.size = 1) + theme(axis.text = element_text(size = 12)) + NoLegend()
dev.off()

pdf("Adult SMG UMAP annotated.pdf", useDingbats = F, width = 3.5, height = 3)
DimPlot(adultsmg, group.by = "CellType", pt.size = 0.5, label = F, label.size = 5, repel = T, cols = colors) + theme(axis.title = element_blank(), axis.text = element_text(size = 12)) + NoLegend()
dev.off()

saveRDS(p30male, file = "P30 male SMG annotated (split from Integrated).rds")
saveRDS(p30female, file = "P30 female SMG annotated (split from Integrated).rds")
saveRDS(adultsmg, file = "Adult SMG annotated (split from Integrated).rds")

# rm(list=ls()) #Clean environment to save space

#################################################################################################
########################## Integration analysis #################################################
#################################################################################################

# If necessary, load previously annotated files since we'll need the metadata
# e12smg <- readRDS("./E12 SMG Annotated (SEURAT v3).rds")
# e14smg <- readRDS("./E14 SMG Annotated (SEURAT v3).rds")
# e16smg <- readRDS("./E16 SMG Annotated (SEURAT v3).rds")
# p1smg <- readRDS("../10X revisions/P1 SMG Annotated (SEURAT v3).rds")
# p30male <- readRDS("../10X revisions/P30 male SMG annotated (split from Integrated).rds")
# p30female <- readRDS("../10X revisions/P30 female SMG annotated (split from Integrated).rds")
# adultsmg <- readRDS("../10X revisions/Adult SMG annotated (split from Integrated).rds")

### Because Integration works in datasets with similar cellular composition (similar cell types are expected)
### we perform integration for embryonic sets and postnatal sets separately.
### Preliminary analysis integrating all datasets resulted in poor resolution and nonsensical results.
### The code for that analysis is provided as reference at the end of this script but it was not used in this manuscript.

e12smg[["sample"]]<-"E12"
e14smg[["sample"]]<-"E14"
e16smg[["sample"]]<-"E16"
p1smg[["sample"]]<-"P1"
p30male[["sample"]]<-"P30M"
p30female[["sample"]]<-"P30F"
adultsmg[["sample"]]<-"Adult"


### Integration of embryonic stages
smg.e.integrated.list <- list(e12smg, e14smg, e16smg)
remove(e12smg, e14smg, e16smg) # remove to save memory

for (i in 1:length(smg.e.integrated.list)) {
  DefaultAssay(smg.e.integrated.list[[i]]) <- "RNA"
  smg.e.integrated.list[[i]] <- NormalizeData(object = smg.e.integrated.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  smg.e.integrated.list[[i]] <- FindVariableFeatures(smg.e.integrated.list[[i]], selection.method = "vst", nfeatures = 2000)
}

integrated.anchors <- FindIntegrationAnchors(smg.e.integrated.list, dims = 1:30)
smg.e.integrated <- IntegrateData(integrated.anchors, dims = 1:30)

remove(smg.e.integrated.list, integrated.anchors) #remove to save memory

# Run the standard workflow for visualization and clustering
smg.e.integrated <- ScaleData(smg.e.integrated, verbose = FALSE)
smg.e.integrated <- RunPCA(smg.e.integrated, npcs = 30)

smg.e.integrated <- FindNeighbors(smg.e.integrated, reduction = "pca", dims = 1:30)
smg.e.integrated <- FindClusters(smg.e.integrated, resolution = 0.6) 
smg.e.integrated <- RunUMAP(smg.e.integrated, reduction = "pca", dims = 1:30)

pdf(file = "Embryonic stages integrated_unsupervised.pdf", useDingbats = F, width = 4.5, height = 4)
DimPlot(smg.e.integrated, reduction = "umap", label = T, group.by = "integrated_snn_res.0.6", label.size = 5) + theme(axis.title = element_blank(), axis.text = element_text(size=12)) + NoLegend()
dev.off()

# make sure identities are sorted as before for consistency
smg.e.integrated <- SetIdent(smg.e.integrated, value = "CellType")
Idents(smg.e.integrated) <- factor(Idents(smg.e.integrated), levels = sort(levels(Idents(smg.e.integrated)),decreasing = F))
smg.e.integrated[["CellType"]]<-Idents(smg.e.integrated) 

stage.colors <- c(
  "#c34d54",
  "#d7cd9e",
  "#019f84",
  "#92d8d4",
  "#0076af",
  "#6b0081"
)
names(stage.colors) <- c("E12", "E14", "E16", "P1", "P30", "Adult")

pdf(file = "Embryonic stages integrated_annotated.pdf", useDingbats = F, width = 6, height = 4.5)
DimPlot(smg.e.integrated, reduction = "umap", group.by = "CellType", label = F, label.size = 4, repel = T, cols = colors) + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
dev.off()

pdf(file = "Embryonic stages integrated_annotated by stage.pdf", useDingbats = F, width = 4.5, height = 4)
DimPlot(smg.e.integrated, reduction = "umap", group.by = "stage", label = F, cols = stage.colors) + theme(axis.title = element_blank(), axis.text = element_text(size=12))
dev.off()

pdf(file = "Embryonic stages integrated_split by stage.pdf", useDingbats = F, width = 3, height = 3)
DimPlot(smg.e.integrated, reduction = "umap", group.by = "stage", label = F, cols = stage.colors, split.by = "stage", ncol = 2) + NoAxes() + NoLegend() 
dev.off()

smg.e.integrated[["celltype.stage"]] <- paste0(smg.e.integrated$stage, "_", smg.e.integrated$CellType)
smg.e.integratedcounts <- as.data.frame(table(smg.e.integrated$celltype.stage))
write.csv(smg.e.integratedcounts, file = "Embryonic cell counts.csv")


saveRDS(smg.e.integrated, file = "Embryonic SMG Integrated (E12-E16).rds")

##########################################
### Integration of postnatal stages 

smg.p.integrated.list <- list(p1smg, p30male, p30female, adultsmg)
remove(p1smg, p30male, p30female, adultsmg) # remove to save memory

for (i in 1:length(smg.p.integrated.list)) {
  DefaultAssay(smg.p.integrated.list[[i]]) <- "RNA"
  smg.p.integrated.list[[i]] <- NormalizeData(object = smg.p.integrated.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  smg.p.integrated.list[[i]] <- FindVariableFeatures(smg.p.integrated.list[[i]], selection.method = "vst", nfeatures = 2000)
}

integrated.anchors <- FindIntegrationAnchors(smg.p.integrated.list, dims = 1:30)
smg.p.integrated <- IntegrateData(integrated.anchors, dims = 1:30)

remove(smg.p.integrated.list, integrated.anchors) #remove to save memory

# Run the standard workflow for visualization and clustering
smg.p.integrated <- ScaleData(smg.p.integrated, verbose = FALSE)
smg.p.integrated <- RunPCA(smg.p.integrated, npcs = 30)

smg.p.integrated <- FindNeighbors(smg.p.integrated, reduction = "pca", dims = 1:30)
smg.p.integrated <- FindClusters(smg.p.integrated, resolution = 0.6) 
smg.p.integrated <- RunUMAP(smg.p.integrated, reduction = "pca", dims = 1:30)

pdf(file = "Postnatal stages integrated_unsupervised.pdf", useDingbats = F, width = 4.5, height = 4)
DimPlot(smg.p.integrated, reduction = "umap", label = T, group.by = "integrated_snn_res.0.6", label.size = 5) + theme(axis.title = element_blank(), axis.text = element_text(size=12)) + NoLegend()
dev.off()


# make sure identities are sorted as before for consistency
smg.p.integrated <- SetIdent(smg.p.integrated, value = "CellType")
Idents(smg.p.integrated) <- factor(Idents(smg.p.integrated), levels = sort(levels(Idents(smg.p.integrated)),decreasing = F))
smg.p.integrated[["CellType"]]<-Idents(smg.p.integrated) 

pdf(file = "Postnatal stages integrated_annotated.pdf", useDingbats = F, width = 6, height = 4.5)
DimPlot(smg.p.integrated, reduction = "umap", group.by = "CellType", label = F, label.size = 5, repel = T, cols = colors) + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
dev.off()

pdf(file = "Postnatal stages integrated_annotated by stage.pdf", useDingbats = F, width = 4.5, height = 4)
DimPlot(smg.p.integrated, reduction = "umap", group.by = "stage", label = F, cols = stage.colors) + theme(axis.title = element_blank(), axis.text = element_text(size=12))
dev.off()

pdf(file = "Postnatal stages integrated_split by stage.pdf", useDingbats = F, width = 3, height = 3)
DimPlot(smg.p.integrated, reduction = "umap", group.by = "stage", label = F, cols = stage.colors, split.by = "stage", ncol = 2) + NoAxes() + NoLegend() 
dev.off()


smg.p.integrated[["celltype.stage"]] <- paste0(smg.p.integrated$stage, "_", smg.p.integrated$CellType)
smg.p.integratedcounts <- as.data.frame(table(smg.p.integrated$celltype.stage))
write.csv(smg.p.integratedcounts, file = "Postnatal cell counts.csv")

saveRDS(smg.p.integrated, file = "Postnatal SMG Integrated (P1-Adult).rds")

## subset and export epithelial clusters for downstream analysis
e.epithelium <- subset(smg.e.integrated, idents = c("Basal duct", "End bud", "Krt19+ duct", "Myoepithelial"))

p.epithelium <- subset(smg.p.integrated, idents = c("Acinar", "Ascl3+ duct", "Basal duct", "Bpifa2+", 
                                                    "Bpifa2+ Proacinar", "GCT", "Intercalated duct", 
                                                    "Krt19+ duct", "Myoepithelial", "Smgc+", "Smgc+ Proacinar",
                                                    "Striated duct", "Mitotic cells"))

pdf("Epithelial cells highlight.pdf", useDingbats = F, width = 4.5, height = 3)
DimPlot(smg.e.integrated, cells.highlight = WhichCells(object = e.epithelium)) + theme(axis.title = element_blank())
DimPlot(smg.p.integrated, cells.highlight = WhichCells(object = p.epithelium)) + theme(axis.title = element_blank())
dev.off()

saveRDS(e.epithelium, file = "Embryonic epithelium.rds")
saveRDS(p.epithelium, file = "Postnatal epithelium.rds")

### Epcam plots
library(RColorBrewer)

pdf("Epcam Plots.pdf", useDingbats = F, width = 2.2, height = 2.2)
FeaturePlot(e12smg, features="Epcam", cols = c("bisque", "green4"), min.cutoff = "q20", max.cutoff = "q80") + theme(axis.text = element_blank(), axis.title = element_blank(), plot.title = element_blank(), axis.ticks = element_blank()) + NoLegend()
FeaturePlot(e14smg, features="Epcam", cols = c("bisque", "green4"), min.cutoff = "q20", max.cutoff = "q80") + theme(axis.text = element_blank(), axis.title = element_blank(), plot.title = element_blank(),  axis.ticks = element_blank()) + NoLegend()
FeaturePlot(e16smg, features="Epcam", cols = c("bisque", "green4"), min.cutoff = "q20", max.cutoff = "q80") + theme(axis.text = element_blank(), axis.title = element_blank(), plot.title = element_blank(), axis.ticks = element_blank()) + NoLegend()
FeaturePlot(p1smg, features="Epcam", cols = c("bisque", "green4"), min.cutoff = "q20", max.cutoff = "q80") + theme(axis.text = element_blank(), axis.title = element_blank(), plot.title = element_blank(), axis.ticks = element_blank()) + NoLegend()
FeaturePlot(p30subset, features="Epcam", cols = c("bisque", "green4"), min.cutoff = "q20", max.cutoff = "q80") + theme(axis.text = element_blank(), axis.title = element_blank(), plot.title = element_blank(), axis.ticks = element_blank()) + NoLegend()
FeaturePlot(smgP.integrated, features="Epcam", cols = c("bisque", "green4"), min.cutoff = "q20", max.cutoff = "q80") + theme(axis.text = element_blank(), axis.title = element_blank(), plot.title = element_blank(), axis.ticks = element_blank()) + NoLegend()
dev.off()

pdf("Epcam Plots2.pdf", useDingbats = F, width = 4.5, height = 4)
FeaturePlot(smg.e.integrated, features="Epcam", cols = c("bisque", "green4"), min.cutoff = "q20", max.cutoff = "q80") + theme(axis.title = element_blank(), axis.text = element_text(size=12), plot.title = element_text(face = "italic", hjust = 0.05, vjust = -5, size = 14)) + NoLegend()
FeaturePlot(smg.p.integrated, features="Epcam", cols = c("bisque", "green4"), min.cutoff = "q20", max.cutoff = "q80") + theme(axis.title = element_blank(), axis.text = element_text(size=12), plot.title = element_text(face = "italic", hjust = 0.05, vjust = -5, size = 14)) + NoLegend()
dev.off()

pdf("Epcam Plots3.pdf", useDingbats = F, width = 3, height = 2.5)
FeaturePlot(e12smg, features="Epcam", cols = c("bisque", "green4"), min.cutoff = "q20", max.cutoff = "q80") + theme(axis.text = element_blank(), axis.title = element_blank(), plot.title = element_text(face = "italic", vjust = -5, hjust = 0.05, size=8))
dev.off()

## END OF SCRIPT ##


