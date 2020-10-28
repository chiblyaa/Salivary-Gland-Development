# This script will perform analysis associated with Figures 5-7
# First we'll perform integratio of all datasets for visualization
# and trajectory analysis. Then, we perform differential expression
# analysis between specific populations and cross-reference with a 
# database of transcription factors. Lastly, we perform correlation
# analysis for Bhlha15 to identify genes associated with acinar lineage.

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(Hmisc)
library(cowplot)
library(viridis)


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

### Load epithelium with updated celltype annotations based on analysis of epithelium
e12smg <- readRDS("../E12-E16 epithelium/E12 epithelium qcd_annotated.rds")
e14smg <- readRDS("../E12-E16 epithelium/E14 epithelium qcd_annotated.rds")
e16smg <- readRDS("../E12-E16 epithelium/E16 epithelium qcd_annotated.rds")
p1smg <- readRDS("../P1-Adult/p1 epithelium qcd_annotated.rds")
p30smg <- readRDS("../P1-Adult/p30 epithelium qcd_annotated.rds")
adultsmg <- readRDS("../P1-Adult/adult epithelium qcd_annotated.rds")

## Perform Integration for all stages 
smg.all.integrated.list <- list(e12smg, e14smg, e16smg, p1smg, p30smg, adultsmg)

for (i in 1:length(smg.all.integrated.list)) {
  DefaultAssay(smg.all.integrated.list[[i]]) <- "RNA"
  smg.all.integrated.list[[i]] <- NormalizeData(object = smg.all.integrated.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  smg.all.integrated.list[[i]] <- FindVariableFeatures(smg.all.integrated.list[[i]], selection.method = "vst", nfeatures = 2000)
}

smg.all.integrated.anchors <- FindIntegrationAnchors(smg.all.integrated.list, dims = 1:30)
smg.all.integrated <- IntegrateData(smg.all.integrated.anchors, dims = 1:30)

remove(smg.all.integrated.list, smg.all.integrated.anchors)

# Run the standard workflow for visualization and clustering
smg.all.integrated <- ScaleData(smg.all.integrated, verbose = FALSE)
smg.all.integrated <- RunPCA(smg.all.integrated, npcs = 30)

smg.all.integrated <- FindNeighbors(smg.all.integrated, reduction = "pca", dims = 1:30)
smg.all.integrated <- FindClusters(smg.all.integrated, resolution = 0.6)
smg.all.integrated <- RunUMAP(smg.all.integrated, reduction = "pca", dims = 1:30)

DimPlot(smg.all.integrated, reduction = "umap", label = T, group.by = "integrated_snn_res.0.6")
DimPlot(smg.all.integrated, reduction = "umap", label = T, group.by = "celltype.fixed")

#remove remaining non-epithelial cells to avoid artifacts in the trajectory
smg.all.integrated <- readRDS("./SMG all epithelium integrated.rds")

smg.all.integrated <- SetIdent(smg.all.integrated, value = "celltype.fixed")
smg.all.integrated <- subset(smg.all.integrated, idents = c("Mesenchyme", "Undefined"), invert=T)
DimPlot(smg.all.integrated)
smg.all.integrated <- RenameIdents(smg.all.integrated, "Smgc+ Female" = "Smgc+", "Smgc+ Male" = "Smgc+")

#Sort identities for consistency
cell.type.levels <- c("End bud", "Krt19+ duct", "Basal duct", "Myoepithelial", "Bpifa2+ Proacinar", "Smgc+ Proacinar", "Mitotic cells",
                      "Bpifa2+", "Acinar", "Smgc+", "Intercalated duct", "Ascl3+ duct", "GCT", "Striated duct")

Idents(smg.all.integrated) <- factor(Idents(smg.all.integrated), levels = cell.type.levels)
smg.all.integrated[["trajectory.identities"]]<-Idents(smg.all.integrated) 
DimPlot(smg.all.integrated, reduction = "umap", group.by = "trajectory.identities", label = T, label.size = 4, repel = T)

pdf(file = "Integrated E12-Adult epithelium UMAPS.pdf", useDingbats = F, width = 5.5, height = 4.5)
DimPlot(smg.all.integrated, reduction = "umap", group.by = "stage", label = F, cols = stage.colors)
DimPlot(smg.all.integrated, reduction = "umap", group.by = "integrated_snn_res.0.6", label = T) + NoLegend()
dev.off()

pdf(file = "Integrated E12-Adult epithelium UMAPS-2.pdf", useDingbats = F, width = 6.5, height = 4.5)
DimPlot(smg.all.integrated, reduction = "umap", group.by = "trajectory.identities", label = F, cols = colors)
dev.off()

saveRDS(smg.all.integrated, file = "SMG all epithelium integrated.rds")


### TRAJECTORY ANALYSIS ###
library(dyno)
library(tidyverse)
library(ggraph)
options(future.globals.maxSize = 2000 * 1024^2)

DimPlot(smg.all.integrated)


epithelium <- smg.all.integrated
DimPlot(epithelium, group.by = "trajectory.identities")
epithelium[["celltype.stage"]] <- paste0(epithelium@meta.data$stage, "_", epithelium@meta.data$trajectory.identities)
epithelium <- SetIdent(epithelium, value = "celltype.stage")
DimPlot(epithelium)

#### create dataset matrix for dynverse package
object_counts <- Matrix::t(as(as.matrix(epithelium@assays$RNA@data), 'sparseMatrix')) ## Note: Integrated SEURAT object does not store "counts" as these are from each individual file prior to integration
object_expression <- Matrix::t(as(as.matrix(epithelium@assays$integrated@data), 'sparseMatrix')) #Use integrated expression values 
dataset<- wrap_expression(
  counts = object_counts, 
  expression = object_expression
)

guidelines <- guidelines_shiny(dataset)
# Reproduces the guidelines as created in the shiny app
answers <- dynguidelines::answer_questions(
  multiple_disconnected = FALSE, 
  expect_topology = FALSE, 
  expected_topology = NULL, 
  expect_cycles = FALSE, 
  expect_complex_tree = TRUE, 
  n_cells = 9777, 
  n_features = 2000, ### only 2000 features are preserved in Integrated data 
  memory = "100GB", 
  prior_information = c("start_id", "end_id", "features_id", "dimred"), 
  docker = TRUE
)
guidelines <- dynguidelines::guidelines(answers = answers)

### create annotation files

epithelium <- SetIdent(epithelium, value = "stage")
e12IDs <- WhichCells(subset(epithelium, idents = "E12"))
endIDs <- WhichCells(subset(epithelium, idents = "Adult"))

groupsids <- as.data.frame(epithelium@meta.data$trajectory.identities)
groupsids$cellID <- colnames(epithelium)
colnames(groupsids) <- c("group_id", "cell_id")

stageids <- as.data.frame(epithelium@meta.data$stage)
stageids$cellID <- colnames(epithelium)
colnames(stageids) <- c("group_id", "cell_id")

clusterids <- as.data.frame(epithelium@meta.data$integrated_snn_res.0.6)
clusterids$cellID <- colnames(epithelium)
colnames(clusterids) <- c("group_id", "cell_id")

cellstage <- as.data.frame(epithelium@meta.data$celltype.stage)
cellstage$cellID <- colnames(epithelium)
colnames(cellstage) <- c("group_id", "cell_id")
cellstage$group_id <- as.character(cellstage$group_id)
cellstage$cell_id <- as.character(cellstage$cell_id)


### calculate trajectory
methods_selected <- guidelines$methods_selected 
dataset <- add_prior_information(dataset = dataset, start_id = e12IDs, end_id = endIDs, start_n = 1, groups_id = groupsids)
dataset <- add_grouping(dataset, epithelium@meta.data$trajectory.identities)

model <- infer_trajectory(dataset, method =  first(methods_selected), give_priors = c("start_id", "end_id"), seed = 15)

pdf("Pseudotime plots.pdf", useDingbats = F, width = 5,height = 5)
plot_dimred(model, grouping = groupsids, color_density = "grouping", color_trajectory = "nearest", color_cells = "grouping", label_milestones = F, size_cells = 0.4, density_cutoff_label = 1, size_trajectory = 1,bw = 0.15, border_radius_percentage = 0.5) + scale_color_manual(values = colors, aesthetics = "colour") + scale_fill_manual(values = colors, aesthetics = "fill") +NoLegend()
plot_dimred(trajectory = model,size_cells = 1, color_cells = "pseudotime", pseudotime = calculate_pseudotime(model), color_density = "grouping", grouping = stageids, label_milestones = F, density_cutoff = 0.1) +scale_discrete_manual(aesthetics = "fill", values = stage.colors) 
dev.off()

model <- model %>% add_root(root_milestone_id = "5") ## milestone 5 was predicted as the start point and it also containes all the cells from E12

pdf("trajectory tree plots.pdf", useDingbats = F, width = 15, height = 7)
#colored by cell type
plot_dendro(trajectory = model,  size_cells = 3, diag_offset = 0.08, grouping = groupsids, color_cells = "grouping", y_offset = 0.6, color_milestones = "grouping", milestones = model$milestone_ids,border_radius_percentage = 0.1) + geom_edge_link() + scale_color_manual(values = colors, aesthetics = "colour") + scale_fill_manual(values = colors, aesthetics = "fill")+ NoLegend()
#colored by stage
plot_dendro(trajectory = model, size_cells = 3, diag_offset = 0.08, alpha_cells = 1, grouping = stageids, color_cells = "grouping", y_offset = 0.6, color_milestones = "none", milestones = model$milestone_ids,border_radius_percentage = 0.1) + geom_edge_link() + scale_color_manual(values = stage.colors, aesthetics = "colour") + scale_fill_manual(values = stage.colors, aesthetics = "fill") + theme(legend.position = "right")
#colored by pseudotime
plot_dendro(trajectory = model, size_cells = 3, diag_offset = 0.08, alpha_cells = 1,  pseudotime = calculate_pseudotime(model), color_cells = "pseudotime", y_offset = 0.6, color_milestones = "none", milestones = model$milestone_ids,border_radius_percentage = 0.1) + geom_edge_link() + theme(legend.position = "right")
dev.off()

### LINEAGE TRACING MARKERS ###
pdf("Lineage tracing markers - trajectory.pdf", useDingbats = F,width = 14,height = 7.5)
plot_dendro(trajectory = model,feature_oi = "Krt14", expression_source = dataset$expression, size_cells = 0.5, alpha_cells = 1,  y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() + scale_fill_viridis(option = "magma", begin = 0.5, end = 1)+scale_color_viridis(option = "magma",begin = 0.5, end = 1, direction = -1) +ggtitle("Krt14")
plot_dendro(trajectory = model,feature_oi = "Acta2", expression_source = dataset$expression, size_cells = 2, alpha_cells = 0.8,  y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() + scale_fill_viridis(option = "magma", begin = 0.5, end = 1)+scale_color_viridis(option = "magma",begin = 0.5, end = 1, direction = -1) +ggtitle("Acta2")
plot_dendro(trajectory = model,feature_oi = "Sox10", expression_source = dataset$expression, size_cells = 2, alpha_cells = 0.8,  y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() + scale_fill_viridis(option = "magma", begin = 0.5, end = 1)+scale_color_viridis(option = "magma",begin = 0.5, end = 1, direction = -1) +ggtitle("Sox2")
plot_dendro(trajectory = model,feature_oi = "Kit", expression_source = dataset$expression, size_cells = 2, alpha_cells = 0.8,  y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() + scale_fill_viridis(option = "magma", begin = 0.5, end = 1)+scale_color_viridis(option = "magma",begin = 0.5, end = 1, direction = -1) +ggtitle("Kit")
plot_dendro(trajectory = model,feature_oi = "Bhlha15", expression_source = dataset$expression, size_cells = 2, alpha_cells = 0.8,  y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() + scale_fill_viridis(option = "magma", begin = 0.5, end = 1)+scale_color_viridis(option = "magma",begin = 0.5, end = 1, direction = -1) +ggtitle("Bhlha15")
plot_dendro(trajectory = model,feature_oi = "Ascl3", expression_source = dataset$expression, size_cells = 2, alpha_cells = 0.8,  y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() + scale_fill_viridis(option = "magma", begin = 0.5, end = 1)+scale_color_viridis(option = "magma",begin = 0.5, end = 1, direction = -1) +ggtitle("Ascl3")
plot_dendro(trajectory = model,feature_oi = "Krt5", expression_source = dataset$expression, size_cells = 2, alpha_cells = 0.8,  y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() + scale_fill_viridis(option = "magma", begin = 0.5, end = 1)+scale_color_viridis(option = "magma",begin = 0.5, end = 1, direction = -1) +ggtitle("Krt5")
dev.off()

saveRDS(model, file = "Trajectory analysis with Integrated SMG epithelium.rds") # import as "model" to use above


########################################################
####### Analysis of developmental transitions ##########
DefaultAssay(epithelium) <-"RNA"

## download list of transcription factors from Schmeier et al (2016): https://tools.sschmeier.com/tcof/home/
## Load file from transcription factors
Transcription.factors <- read.csv("~/10X Manuscript Projects/10X revisions/Transition Analysis/TranscriptionFactors.csv")
names(Transcription.factors)[1] <- "Gene"

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

### create additional metadata for analyis of developmental transitions
epithelium[["transition.labels"]] <- epithelium$celltype.stage
DimPlot(epithelium, group.by = "transition.labels")
epithelium <- SetIdent(epithelium, value = "transition.labels")

#################### DUCT POPULATIONS - set appropriate identity for comparisons
# extract only populations related to duct lineage for visualization
epithelium <- SetIdent(epithelium, value = "celltype.fixed")
duct.idents <- c("Krt19+ duct", "Basal duct", "Ascl3+ duct",
                 "Striated duct", "GCT", "Intercalated duct")
duct.epithelium <- subset(epithelium, idents = duct.idents)


epithelium[["duct.labels"]] <- "NA"
epithelium[["duct.labels"]] <- Idents(duct.epithelium)
pdf("duct populations.pdf", width = 6.5,height = 5,useDingbats = F)
DimPlot(epithelium, group.by = "duct.labels", cols = colors) + theme(axis.text = element_text(size = 12))
dev.off()

duct.epithelium <- SetIdent(duct.epithelium, value = "transition.labels")

#create a table with the desired comparisons:
duct.pops <- data.frame('1'=c("E12_Krt19+ duct", "E14_Krt19+ duct"),
                          '2'=c("E14_Krt19+ duct", "E16_Krt19+ duct"),
                          '3'=c("E14_Basal duct", "E16_Basal duct"),
                          '4'=c("E16_Krt19+ duct", "P1_Krt19+ duct"),
                          '5'=c("E16_Basal duct", "P1_Basal duct"),
                          '6'=c("P1_Basal duct", "P30_Basal duct"), #5 and 2
                          '7'=c("P1_Krt19+ duct", "P30_Intercalated duct"),
                          '8'=c("P1_Krt19+ duct", "P30_Ascl3+ duct"),
                          '9' = c("P1_Krt19+ duct", "P30_Striated duct"),
                          '10' = c("P1_Krt19+ duct", "P30_GCT"), stringsAsFactors = F)

duct.tests <- as.data.frame(t(as.matrix(duct.pops, rownames.force = T)))
filenames <- c("E12_Krt19 duct to E14_Krt19  duct",
               "E14_Krt19 duct to E16_Krt19  duct",
               "E14_Basal duct to E16_Basal duct",
               "E16_Krt19  duct to P1_Krt19  duct",
               "E16_Basal duct to P1_Basal duct",
               "P1_Basal duct to P30_Basal duct", #5 and 2
               "P1_Krt19  duct to P30_Intercalated duct",
               "P1_Krt19  duct to P30_Ascl3  duct",
               "P1_Krt19  duct to P30_Striated duct",
               "P1_Krt19  duct to P30_GCT")
rownames(duct.tests) <- filenames
names(duct.tests) <- c("Ident1", "Ident2")
duct.tests$Ident1 <- as.character(duct.tests$Ident1)
duct.tests$Ident2 <- as.character(duct.tests$Ident2)

epithelium<-SetIdent(epithelium, value = "transition.labels")
duct.comparisons <- list()
for (j in 1:NROW(duct.tests)) {
  duct.comparisons[[j]] <- transition.analysis(seuratobj = epithelium, ident1 = duct.tests$Ident1[j], ident2 = duct.tests$Ident2[j], 
                                                 tfs = Transcription.factors, 
                                                 cols = stage.colors,
                                                 filename = rownames(duct.tests)[j], up = 5, down = 5)
}

#additional analyses from basal duct to GCT, SD, and acinar
bdacinar <- transition.analysis(seuratobj = epithelium, ident1 = c("P30_Basal duct", "Adult_Basal duct"), ident2 = c("P30_Acinar", "Adult_Acinar"), 
                                tfs = Transcription.factors, 
                                cols = colors,
                                filename = "Basal duct to Acinar", up = 5, down = 5, group.by = "celltype.fixed")
bdsd <- transition.analysis(seuratobj = epithelium, ident1 = c("P30_Basal duct", "Adult_Basal duct"), ident2 = c("P30_Striated duct", "Adult_Striated duct"), 
                                tfs = Transcription.factors, 
                                cols = colors,
                                filename = "Basal duct to Striated duct", up = 5, down = 5, group.by = "celltype.fixed")
bdgct <- transition.analysis(seuratobj = epithelium, ident1 = c("P30_Basal duct", "Adult_Basal duct"), ident2 = c("P30_GCT"), 
                            tfs = Transcription.factors, 
                            cols = colors,
                            filename = "Basal duct to GCT", up = 5, down = 5, group.by = "celltype.fixed")



############## ACINAR LINEAGE POPULATIONS 
# extract only populations related to acinar lineage
epithelium <- SetIdent(epithelium, value = "celltype.fixed")
as.data.frame(table(Idents(epithelium)))
acinar.idents <- c("End bud", "Smgc+ Proacinar", "Smgc+ Male", "Smgc+ Female", "Smgc+",
                   "Bpifa2+ Proacinar", "Acinar", "Bpifa2+")
acinar.epithelium <- subset(epithelium, idents = acinar.idents)
acinar.epithelium <- RenameIdents(acinar.epithelium, "Smgc+ Female" ="Smgc+", "Smgc+ Male" = "Smgc+") ## we do not need to differentiate between male and female for this, we combine for simplicity

## visualization
epithelium[["acinar.labels"]] <- "NA"
epithelium[["acinar.labels"]] <- Idents(acinar.epithelium)
pdf("Acinar populations.pdf", width = 6.5,height = 5,useDingbats = F)
DimPlot(epithelium, group.by = "acinar.labels", cols = colors) + theme(axis.text = element_text(size = 12))
dev.off()

#create a table with the desired comparisons:
acinar.pops <- data.frame('1'=c("E12_End bud", "E14_End bud"),
                          '2'=c("E14_End bud", "E16_End bud"),
                          '3'=c("E16_End bud", "P1_Smgc+ Proacinar"),
                          '4'=c("E16_End bud", "P1_Bpifa2+ Proacinar"),
                          '5'=c("P1_Smgc+ Proacinar", "P30_Smgc+"),
                          '7'=c("P1_Smgc+ Proacinar", "P30_Acinar"),
                          '8'=c("P1_Bpifa2+ Proacinar", "P30_Acinar"), stringsAsFactors = F)

acinar.tests <- as.data.frame(t(as.matrix(acinar.pops, rownames.force = T)))
filenames <- c("E12 to E14 End bud", "E14 to E16 End bud", "E16 End bud to Smgc proacinar",
               "E16 End bud to Bpifa2 proacinar", "Smgc proacinar to Smgc", "Smgc proacinar to P30 acinar",
               "Bpifa2 proacinar to P30 Acinar")
rownames(acinar.tests) <- filenames
names(acinar.tests) <- c("Ident1", "Ident2")
acinar.tests$Ident1 <- as.character(acinar.tests$Ident1)
acinar.tests$Ident2 <- as.character(acinar.tests$Ident2)

epithelium <- SetIdent(epithelium, value = "transition.labels")
DimPlot(epithelium)

acinar.comparisons <- list()
for (j in 1:NROW(acinar.tests)) {
  acinar.comparisons[[j]] <- transition.analysis(seuratobj = epithelium, ident1 = acinar.tests$Ident1[j], ident2 = acinar.tests$Ident2[j], 
                                                 tfs = Transcription.factors, 
                                                 cols = stage.colors,
                                                 filename = rownames(acinar.tests)[j], up = 5, down = 5)
}


################# CORRELATION ANALYSIS FOR ACINAR MARKERS #########################
library(corrplot)
library(RColorBrewer)

# extract expression values and perform correlation analysis for Mist1
DefaultAssay(acinar.epithelium) <- "RNA"
matrix<-GetAssayData(acinar.epithelium)
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Bhlha15",])
correlations<-apply(matrix_mod,1,function(x){cor.test(gene,x)})

df <- data.frame(matrix(unlist(correlations), NROW(acinar.epithelium),  byrow=T),stringsAsFactors=FALSE) # 23774 rows in RNA assay
rownames(df) <- names(correlations)
colnames(df) <- names(correlations$Xkr4)
names(df)[3] <-"pvalue"
names(df)[4] <-"R"
df$pvalue <- as.numeric(df$pvalue)
df$R <- as.numeric(df$R)

write.csv(df, file = "Mist1-correlated genes.csv")
dftfs <- df[rownames(df) %in% Transcription.factors$Gene, ]
write.csv(dftfs, file = "Mist-1 correlated TFs.csv")

## We'll do this in a different way for visualization of the top correlated markers
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
correlations2 <- as.data.frame(correlations)
correlations2$gene <- rownames(correlations2)
correlations3 <- correlations2[correlations2$gene %in% Transcription.factors$Gene, ]
correlations3 <- rbind(correlations3 %>% top_n(10,correlations), correlations3 %>% top_n(-10, correlations))
correlations3 <- correlations3[order(correlations3$correlations), ]

#create matrix subset with transcription factors only to speed up analysis and get correlations between them 
sub.matrix <- subset.matrix(matrix_mod, rownames(matrix_mod) %in% Transcription.factors$Gene)
sub.matrix <- as.matrix(t(sub.matrix))
test.run <- cor(sub.matrix)
test.run.subset <- subset.matrix(test.run, rownames(test.run) %in% correlations3$gene)
test.run.subset <- subset.matrix(test.run.subset, select = rownames(test.run) %in% correlations3$gene)

# custom function to get p-values of correlations between top transcription factors (use for visualization)
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.matrix <- subset.matrix(matrix_mod, rownames(matrix_mod) %in% correlations3$gene)
p.matrix <- as.matrix(t(p.matrix))
p.mat <- cor.mtest(p.matrix)
head(p.mat[, 1:5])

corrplot(test.run.subset, method="circle",tl.col = "black", order = "hclust",
         col = rev(brewer.pal(n=8, name="Spectral")),type = "upper",sig.level = 0.01, p.mat = p.mat, 
         insig = "blank") + theme(panel.background = element_blank())

pdf(file = "Correlation plot for Mist1-correlated transcription factors.pdf", useDingbats = F, width = 5,height = 5)
corrplot(test.run.subset, method="color",bg = "transparent",diag = T,outline = "black",addgrid.col = "black",
         tl.col = "black", order = "FPC", 
         col = rev(brewer.pal(n=8, name="Spectral")),type = "upper",sig.level = 0.01, insig="pch", p.mat = p.mat) 
dev.off()

#another visualization tool
library(heatmaply)
heatmaply_cor(test.run.subset)

acinar.epithelium <- SetIdent(acinar.epithelium, value = "stage")
#Heatmap visualization across dev stages
DefaultAssay(acinar.epithelium) <- "RNA"

library(superheat)
library(viridis)
acinar.gene.matrix <- as.matrix(GetAssayData(acinar.epithelium),slot = "counts") 
acinar.gene.matrix <- subset(acinar.gene.matrix, subset = rownames(acinar.gene.matrix) %in% correlations3$gene)
acinar.epithelium$stage <- factor(acinar.epithelium$stage, levels = c("E12", "E14", "E16", "P1", "P30", "Adult")) # correct order of dev stages

pdf(file = "Heatmap of mist1-correlated top tfs.pdf", useDingbats = F, width = 5, height = 5)
superheat(t(acinar.gene.matrix),
          scale = T,legend = T, pretty.order.cols = T, 
          bottom.label.text.angle = 90, bottom.label.size = 0.5,  
          bottom.label.text.size = 4, left.label.text.size = 0, 
          left.label.col = stage.colors, bottom.label.col = "white",
          grid.hline.size = 0.2, grid.vline.size = 0.2, 
          left.label.size = 0.05, bottom.label.text.alignment = "right", membership.rows = acinar.epithelium$stage,
          heat.pal = rev(brewer.pal(8,name = "Spectral")), force.grid.hline = T, legend.height = 0.2,legend.width = 1, pretty.order.rows = T)
dev.off()

saveRDS(epithelium, "epithelium backup.rds")
