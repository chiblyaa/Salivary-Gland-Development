library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
##SEURAT V3#

#####################################################################################################
###################### ADULT SMG -- LOAD AND DATA AND CLUSTER WITH SEURAT ###########################
#####################################################################################################

#### Load data matrix from P30 ICR and Adult C3H mouse

SMG.ICR.data <- Read10X(data.dir = "~/Desktop/Adult SMG (ICR CH3)/ICRaggr_v3/filtered_gene_bc_matrices_mex/mm10plus/")
SMG.CH3.data <- Read10X(data.dir = "~/Desktop/Adult SMG (ICR CH3)/SMGtest/filtered_feature_bc_matrix/")
ICRsmg <- CreateSeuratObject(SMG.ICR.data, min.cells = 5, min.features = 200)
CH3smg <- CreateSeuratObject(SMG.CH3.data, min.cells = 5, min.features = 200)

##### QC DATA
ICRsmg[["percent.mt"]] <- PercentageFeatureSet(object = ICRsmg, pattern = "^mt-")
CH3smg[["percent.mt"]] <- PercentageFeatureSet(object = CH3smg, pattern = "^mt-")
VlnPlot(object = ICRsmg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = CH3smg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ICRsmg <- subset(x = ICRsmg, subset = nFeature_RNA > 200 &  percent.mt < 10)
CH3smg <- subset(x = CH3smg, subset = nFeature_RNA > 200 &  percent.mt < 10)

### Normalize and add metadata
ICRsmg <- NormalizeData(object = ICRsmg, normalization.method = "LogNormalize", scale.factor = 10000)
CH3smg <- NormalizeData(object = CH3smg, normalization.method = "LogNormalize", scale.factor = 10000)
ICRsmg$strain <- "ICR"
CH3smg$strain <- "CH3"
ICRsmg <- FindVariableFeatures(object = ICRsmg, selection.method = "vst", nfeatures = 2000)
CH3smg <- FindVariableFeatures(object = CH3smg, selection.method = "vst", nfeatures = 2000)

#### merge datasets

adult.smg <- merge(x = ICRsmg, y = CH3smg, add.cell.ids = c("ICR", "C3H"), merge.data = F)
adult.smg <- NormalizeData(adult.smg)
adult.smg <- FindVariableFeatures(adult.smg, selection.method = "vst",nfeatures = 2000)
# Run the standard workflow for visualization and clustering
adult.smg <- ScaleData(object = adult.smg, verbose = T)
adult.smg <- RunPCA(object = adult.smg, npcs = 30, verbose = FALSE)
ElbowPlot(adult.smg)
# t-SNE and Clustering

adult.smg <- FindNeighbors(object = adult.smg, reduction = "pca", dims = 1:12, force.recalc = T)
adult.smg <- FindClusters(adult.smg, resolution = 0.6)
adult.smg <- RunTSNE(object = adult.smg, reduction = "pca", dims = 1:12)

DimPlot(object = adult.smg, reduction = "tsne", group.by = "strain")
DimPlot(object = adult.smg, reduction = "tsne", label = TRUE)

### Identify cell cluster markers
allcellmarkers <- FindAllMarkers(adult.smg, logfc.threshold = 0.25, only.pos = T)
allcellmarkers <- allcellmarkers[allcellmarkers$p_val_adj<0.05,]
topmarkers <- allcellmarkers %>% group_by(cluster) %>% top_n(5,avg_logFC)
write.table(x = allcellmarkers, file = "Adult SMG Unsupervised cluster markers.txt", sep = "\t", row.names = F)
DotPlot(adult.smg, features = rev(unique(topmarkers$gene)), dot.scale = 6, group.by = "seurat_clusters", cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))

### Generate plot with unsupervised clusters and look at strain integration
DimPlot(adult.smg, pt.size = 1,reduction = "tsne", label.size = 6,repel = F, label = T)
DimPlot(adult.smg, pt.size = 1,reduction = "tsne", group.by = "strain")
FeaturePlot(adult.smg, features = "Egf", pt.size = 1, min.cutoff = "q10")

### Annotate clusters with cell ID and save under metadata "CellType"
adult.smg <- RenameIdents(object = adult.smg, `0` = "Acinar", `1` = "Smgc+", `2` = "Bpifa2+ Acinar", `3` = "Endothelial", `4` = "Smgc+", `5` = "Striated duct", `6` = "Acinar", `7` = "Basal duct", `8` = "Intercalated duct", `9` = "GCT", `10` = "Immune", `11` = "GCT", `12` = "Immune", `13` = "Acinar", `14` = "Myoepithelial", '15' = "Ascl3", '16'="Immune", '17'="Erythroid")
adult.smg[["CellType"]] <- Idents(adult.smg)
adult.smg[["stage"]] <- "Adult"

### Generate TNSE for Figure and export file for downstream analysis
DimPlot(adult.smg,group.by = "CellType", pt.size = 1,reduction = "tsne", label.size = 6,repel = F, label = T)

adult.smg<- SetIdent(adult.smg, value = "CellType")
adult.smg.celltype.markers <- FindAllMarkers(adult.smg, only.pos = T, logfc.threshold = 0.25)
adult.smg.celltype.markers <-adult.smg.celltype.markers[adult.smg.celltype.markers$p_val_adj <0.05, ]
write.table(adult.smg.celltype.markers, file = "Adult SMG cell type markers.txt", sep = "\t", row.names = F)

adult.cellmarkers.top5 <- adult.smg.celltype.markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
DotPlot(adult.smg, features = rev(unique(adult.cellmarkers.top5$gene)), dot.scale = 6, group.by = "CellType", cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))


saveRDS(adult.smg, file = "Adult SMG (ICR-CH3) Annotated (SEURAT v3).rds")


### integrated analysis - for comparison

all.features <- unique(c(rownames(ICRsmg), rownames(CH3smg)))
smg.anchors <- FindIntegrationAnchors(object.list = c(ICRsmg, CH3smg), dims = 1:30)
smg.integrated <- IntegrateData(anchorset = smg.anchors, dims = 1:30, features.to.integrate = all.features)
DefaultAssay(pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
smg.integrated <- RunPCA(smg.integrated, npcs = 30, verbose = FALSE)
smg.integrated <- RunUMAP(smg.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(smg.integrated, reduction = "umap", group.by = "strain")
p2 <- DimPlot(smg.integrated, reduction = "umap", group.by = "CellType", label = TRUE, 
              repel = TRUE) + NoLegend()
plot_grid(p1, p2)

###label clusters based on previous markers
### code not shown because it was performed with previous version of seurat no longer compatible with object

DimPlot(smg.integrated)
smg.integrated <- RenameIdents(smg.integrated, 'Immune 1' = "Immune", 'Immune 2' = "Immune")
smg.integrated[["Epi.CellType"]] <- Idents(smg.integrated)

smg.integrated <- SetIdent(smg.integrated,value = "Epi.CellType")
smgIntegrated.markers <- FindAllMarkers(smg.integrated, logfc.threshold = 0.25, only.pos = T)
adult.smg.celltype.markers <-adult.smg.celltype.markers[adult.smg.celltype.markers$p_val_adj <0.05, ]
write.table(adult.smg.celltype.markers, file = "Adult SMG cell type markers.txt", sep = "\t", row.names = F)

adult.cellmarkers.top5 <- adult.smg.celltype.markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
DotPlot(adult.smg, features = rev(unique(adult.cellmarkers.top5$gene)), dot.scale = 6, group.by = "CellType", cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))

smg.integrated <- SetIdent(smg.integrated, value = "seurat_clusters")
smgIntegrated.markers <- FindAllMarkers(smg.integrated, logfc.threshold = 0.25, only.pos = T)

saveRDS(smg.integrated, file = "Adult SMG (Integrated ICR-CH3) Annotated (SEURAT v2).rds")



