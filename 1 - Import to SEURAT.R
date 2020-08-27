library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(Hmisc)
library(circlize)

e12smg.data <- Read10X(data.dir = "~/CellRanger Files/E12_V3/outs/filtered_gene_bc_matrices/mm10plus/")
e14smg.data <- Read10X(data.dir = "~/CellRanger Files/E14_V3/outs/filtered_gene_bc_matrices/mm10plus/")
e16smg.data <- Read10X(data.dir = "~/CellRanger Files/SalGland_E16_v2/filtered_gene_bc_matrices/mm10plus/")
p1smg.data <- Read10X(data.dir = "~/CellRanger Files/SalGland_P1_v3/filtered_gene_bc_matrices/mm10plus/")

#### THE FOLLOWING CODE IS PROVIDED AS REFERENCE FOR HOW TO ANALYZE DATA USING SEURAT V3 DEFAULT PIPELINE
#### ORIGINAL ANALYSIS WAS PERFORMED USING SEURAT V2, WHICH IS NOW OUTDATED.

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
DimPlot(object = e12smg, reduction = "tsne", pt.size = 1,label = T,label.size = 6)

#### Identify cluster markers for annotation

e12markers <- FindAllMarkers(e12smg, only.pos = T,logfc.threshold = 0.25)
e12markers.top5 <- e12markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.table(e12markers, file = "E12 SMG Unsupervised markers (SEURAT).txt", sep = "\t", row.names = F)
DotPlot(e12smg, features = unique(e12markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

e12smg[["stage"]] <- "E12"
e12smg <- RenameIdents(e12smg, '0' = "Mesenchyme",'1' = "Mesenchyme",'2' = "Mesenchyme",'3' = "Mesenchyme",'4' = "Mesenchyme",'5' = "Mesenchyme",'6' = "Mesenchyme", 
                       '7'= "End bud", '8' = "Krt19+ duct", '9'= "Nerves", '10'="Smooth muscle", '11' = "Erythroid", '12'="Endothelial", '13'="Immune")
e12smg[["CellType"]] <- Idents(e12smg)
DimPlot(e12smg, group.by = "CellType", label = T, label.size = 6, pt.size = 1)

e12smg <- SetIdent(e12smg, value = "CellType")
e12cellcounts <- as.data.frame(table(Idents(e12smg)))
e12.cell.markers <- FindAllMarkers(e12smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
e12.cell.markers.top5 <- e12.cell.markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.table(e12.cell.markers, file = "E12 SMG Cell Type markers (SEURAT).txt", sep = "\t", row.names = F)
DotPlot(e12smg, features = unique(e12.cell.markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "CellType") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

saveRDS(e12smg, file = "E12 SMG Annotated (SEURAT v3).rds")

#############################################################################################
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
DimPlot(object = e14smg, reduction = "tsne", pt.size = 1,label = T,label.size = 6)

e14markers <- FindAllMarkers(e14smg, only.pos = T,logfc.threshold = 0.25)
e14markers.top5 <- e14markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.table(e14markers, file = "e14 SMG Unsupervised markers (SEURAT).txt", sep = "\t", row.names = F)
DotPlot(e14smg, features = unique(e14markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

e14smg[["stage"]] <- "E14"
e14smg <- RenameIdents(e14smg, '0' = "Mesenchyme",'1' = "Mesenchyme",'3' = "Mesenchyme",'14' = "Mesenchyme",'8' = "Mesenchyme",
                       '19' = "Endothelial",'16' = "Nerves",'17'="Immune",'18'="Immune",'11' = "Erythroid",'15' = "Erythroid", '13'="Glial cells",
                       '7'= "Krt19+ duct", '6' = "Krt19+ duct", '5'="End bud", '2'="End bud",'9'="End bud",'4'="End bud",'12'="End bud",'10'="End bud")
e14smg[["CellType"]] <- Idents(e14smg)
DimPlot(e14smg, group.by = "CellType", label = T, label.size = 6, pt.size = 1)

e14smg <- SetIdent(e14smg, value = "CellType")
e14cellcounts <- as.data.frame(table(Idents(e14smg)))
e14.cell.markers <- FindAllMarkers(e14smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
e14.cell.markers.top5 <- e14.cell.markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.table(e14.cell.markers, file = "e14 SMG Cell Type markers (SEURAT).txt", sep = "\t", row.names = F)
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
DimPlot(object = e16smg, reduction = "tsne", pt.size = 1,label = T,label.size = 6)

e16markers <- FindAllMarkers(e16smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
write.table(e16markers, file = "e16 SMG Unsupervised markers (SEURAT).txt", sep = "\t", row.names = F)
e16markers.top5 <- e16markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
DotPlot(e16smg, features = unique(e16markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

e16smg[["stage"]] <- "e16"
e16smg <- RenameIdents(e16smg, '0' = "Mesenchyme",'1' = "Mesenchyme",'2' = "Mesenchyme",'3' = "Mesenchyme",'6' = "Mesenchyme",'8' = "Mesenchyme",'15' = "Mesenchyme",'16' = "Mesenchyme",
                       '12' = "Erythroid", '11'="Immune",'18'="Immune", '21' = "Nerves", '20'="Glial cells", '19' = "Immune", '17'="Smooth muscle", '14'="Basal duct", '13'="Krt19+ duct",
                       '10'= "Endothelial", '9' = "Myoepithelial", '5'="Erythroid", '4'="End bud",'7'="Erythroid")
e16smg[["CellType"]] <- Idents(e16smg)
DimPlot(e16smg, group.by = "CellType", label = T, label.size = 6, pt.size = 1) +NoLegend()

e16smg <- SetIdent(e16smg, value = "CellType")
e16cellcounts <- as.data.frame(table(Idents(e16smg)))
e16.cell.markers <- FindAllMarkers(e16smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
e16.cell.markers.top5 <- e16.cell.markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.table(e16.cell.markers, file = "e16 SMG Cell Type markers (SEURAT).txt", sep = "\t", row.names = F)
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

p1markers <- FindAllMarkers(p1smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
write.table(p1markers, file = "p1 SMG Unsupervised markers (SEURAT).txt", sep = "\t", row.names = F)
p1markers.top5 <- p1markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
DotPlot(p1smg, features = unique(p1markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "seurat_clusters") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

p1smg[["stage"]] <- "p1"
p1smg <- RenameIdents(p1smg, '0' = "Bpifa2+ Proacinar",'9' = "Bpifa2+ Proacinar",'10' = "Bpifa2+ Proacinar",'1' = "Smgc+ Proacinar", '2' = "Smgc+ Proacinar",'3' = "Smgc+ Proacinar",
                      '17'= "Endothelial", '5'="Krt19+ duct", '6' = "Smgc+ Proacinar", '7' = "Basal duct", '8'= "Mitotic cells", '11' = "Bpifa2+ Proacinar", '12'="Krt19+ duct",
                      '15' = "Basal duct", '14' = "Myoepithelial", '16'="Immune", '4' = "Mesenchyme", '13'="Mesenchyme", '18' ="Erythroid")
p1smg[["CellType"]] <- Idents(p1smg)
DimPlot(p1smg, group.by = "CellType", label = T, label.size = 6, pt.size = 1) +NoLegend()

p1smg <- SetIdent(p1smg, value = "CellType")
p1cellcounts <- as.data.frame(table(Idents(p1smg)))
p1.cell.markers <- FindAllMarkers(p1smg, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
p1.cell.markers.top5 <- p1.cell.markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.table(p1.cell.markers, file = "p1 SMG Cell Type markers (SEURAT).txt", sep = "\t", row.names = F)
DotPlot(p1smg, features = unique(p1.cell.markers.top5$gene), cols = "Spectral", dot.scale = 6,group.by = "CellType") +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

saveRDS(p1smg, file = "p1 SMG Annotated (SEURAT v3).rds")

###### END OF SCRIPT ######


