library(dplyr)
library(Seurat)
library(Hmisc)

########### LOAD PREVIOUSLY SAVED FILE WITH EPITHELIAL CLUSTERS FROM ALL STAGES ###########################

e12epithelium <- readRDS("~/1 - Separate epithelial clusters/E12 Epithelium annotated (SEURAT v3).rds")
e14epithelium <- readRDS("~/1 - Separate epithelial clusters/E14 Epithelium annotated (SEURAT v3).rds")
e16epithelium <- readRDS("~/1 - Separate epithelial clusters/E16 Epithelium annotated (SEURAT v3).rds")
p1epithelium <- readRDS("~/1 - Separate epithelial clusters/P1 Epithelium annotated (SEURAT v3).rds")
adult.epithelium <- readRDS("~/1 - Separate epithelial clusters/Adult (C3H ICR) Epithelium annotated (SEURAT v3).rds")

########### Remove remaining non-epithelial cells from each file ###########################

DimPlot(e12epithelium, group.by = "Epi.CellType")
e12epithelium <- SetIdent(e12epithelium, value = "Epi.CellType")
e12epithelium <- subset(e12epithelium,idents = c("End bud", "Krt19+ duct"))
e12epithelium[["Epi.CellType.Stage"]] <- paste0("E12_", e12epithelium$Epi.CellType)
e12epithelium[["Cluster.Stage"]] <- paste0("E12_", e12epithelium$seurat_clusters)

DimPlot(e14epithelium, group.by = "Epi.CellType")
e14epithelium <- SetIdent(e14epithelium, value = "Epi.CellType")
e14epithelium <- subset(e14epithelium, idents = c("End bud", "Duct"))
e14epithelium <- RenameIdents(e14epithelium, 'Duct' = "Krt19+ duct")
e14epithelium$Epi.CellType <- Idents(e14epithelium)
e14epithelium[["Epi.CellType.Stage"]] <- paste0("E14_", e14epithelium$Epi.CellType)
e14epithelium[["Cluster.Stage"]] <- paste0("E14_", e14epithelium$seurat_clusters)

DimPlot(e16epithelium, group.by = "Epi.CellType")
e16epithelium <- SetIdent(e16epithelium, value = "Epi.CellType")
e16epithelium <- subset(e16epithelium,idents = c("End bud", "Krt19+ duct", "Basal duct", "Myoepithelial"))
e16epithelium[["Epi.CellType.Stage"]] <- paste0("E16_", e16epithelium$Epi.CellType)
e16epithelium[["Cluster.Stage"]] <- paste0("E16_", e16epithelium$seurat_clusters)

DimPlot(p1epithelium, group.by = "Epi.CellType")
p1epithelium <- RenameIdents(p1epithelium, "Bpifa2+ Smgc+ Proacinar" = "Bpifa2+ Proacinar", 'Duct cluster 4' = "Krt19+ duct")
p1epithelium[["Epi.CellType"]] <- Idents(p1epithelium)
DimPlot(p1epithelium, label = T, label.size = 6, group.by = "Epi.CellType")
p1epithelium[["Epi.CellType.Stage"]] <- paste0("P1_", p1epithelium$Epi.CellType)
p1epithelium[["Cluster.Stage"]] <- paste0("P1_", p1epithelium$seurat_clusters)

DimPlot(adult.epithelium, label = T, label.size = 6,group.by = "Epi.CellType")
adult.epithelium[["Epi.CellType.Stage"]] <- paste0("Adult_", adult.epithelium$Epi.CellType)
adult.epithelium[["Cluster.Stage"]] <- paste0("Adult_", adult.epithelium$seurat_clusters)

########### Integrate epithelium from all stages into single file ###########################

epithelium <- merge(e12epithelium, y = c(e14epithelium, e16epithelium, p1epithelium, adult.epithelium), add.cell.ids = c("E12", "E14", "E16", "P1", "Adult"), merge.data = T)
remove(e12epithelium, e14epithelium, e16epithelium, p1epithelium, adult.epithelium)

###########################################################################################################
### RE-CLUSTER EPITHELIUM #################################################################################

epithelium <- NormalizeData(epithelium)
epithelium <- FindVariableFeatures(object = epithelium, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = epithelium), 10)
plot1 <- VariableFeaturePlot(object = epithelium)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(x = epithelium)
epithelium <- ScaleData(object = epithelium, features = all.genes)
epithelium <- RunPCA(object = epithelium, features = VariableFeatures(object = epithelium))
ElbowPlot(object = epithelium)

epithelium <- FindNeighbors(object = epithelium, dims = 1:12)
epithelium <- FindClusters(object = epithelium, resolution = 0.5)
epithelium <- RunTSNE(object = epithelium, dims = 1:12)

DimPlot(object = epithelium, reduction = "tsne", group.by = "seurat_clusters", label = T, pt.size = 1, label.size = 8, repel = F) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) + NoLegend()
DimPlot(object = epithelium, reduction = "tsne", group.by = "Epi.CellType", label = F, pt.size = 1, label.size = 5, repel = F, cols = c("navy", "hotpink", "darkred", "slateblue1", "slateblue3", "slateblue4", "firebrick2", "skyblue2", "goldenrod2", "purple", "firebrick", "forestgreen", "mediumpurple2", "pink2", "turquoise", "gray", "gray")) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) 
DimPlot(object = epithelium, reduction = "tsne", group.by = "Epi.CellType.Stage", label = F, pt.size = 1, label.size = 5, repel = F, cols = c("navy", "hotpink", "darkred","slateblue1", "goldenrod2", "purple2", "forestgreen", "pink2", "turquoise2", "lightskyblue1", "orange2", "lightskyblue3", "palevioletred3", "firebrick", "lightskyblue4", "orange", "green", "orangered3", "slateblue2", "darkorange2", "springgreen4", "orchid")) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) 
DimPlot(object = epithelium, reduction = "tsne", group.by = "stage", label = F, label.size = 6, pt.size = 1) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) 


saveRDS(epithelium, file = "Merged epithelium from all stages (Seurat v3).rds")

