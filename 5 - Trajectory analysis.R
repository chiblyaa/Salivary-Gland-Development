library(Seurat)
library(dyno)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggraph)
library(RColorBrewer)

### Load data with annotated epithelium (merged with all stages)
epithelium <- readRDS("~/2 - Integrate epithelial clusters/Merged epithelium from all stages (Seurat v3).rds")
epithelium <- SetIdent(epithelium, value="Epi.CellType.Stage")
DimPlot(epithelium)

#### create dataset matrix for dynverse package
object_counts <- Matrix::t(as(as.matrix(epithelium@assays$RNA@counts), 'sparseMatrix'))
object_expression <- Matrix::t(as(as.matrix(epithelium@assays$RNA@data), 'sparseMatrix'))
dataset<- wrap_expression(
  counts = object_counts, 
  expression = object_expression
)

guidelines <- guidelines_shiny(dataset)

## Reproduces the guidelines as created in the shiny app
answers <- dynguidelines::answer_questions(
  multiple_disconnected = FALSE, 
  expect_topology = TRUE, 
  expected_topology = "tree", 
  n_cells = 9786, 
  n_features = 23737, 
  memory = "100GB", 
  prior_information = c("start_id", "end_id", "end_n", "start_n", "leaves_n", "groups_n", "features_id"), 
  docker = TRUE
)
guidelines <- dynguidelines::guidelines(answers = answers) 

### create annotation files
e12IDs <- dataset$cell_ids[grep("E12" ,dataset$cell_ids)]
AdultIDs <- dataset$cell_ids[grep("Adult" ,dataset$cell_ids)]
groupsids <- as.data.frame(epithelium@meta.data$Epi.CellType)
groupsids$cellID <- colnames(epithelium)
colnames(groupsids) <- c("group_id", "cell_id")
stageids <- as.data.frame(epithelium@meta.data$stage)
stageids$cellID <- colnames(epithelium)
colnames(stageids) <- c("group_id", "cell_id")

######## Trajectory analysis
methods_selected <- guidelines$methods_selected 
dataset <- add_prior_information(dataset, start_id = e12IDs, end_id = AdultIDs, start_n = 1,verbose = T, groups_id = groupsids)
dataset <- add_grouping(dataset, epithelium@meta.data$Epi.CellType.Stage)
model <- infer_trajectory(dataset, first(methods_selected))

model$milestone_network
head(model$progressions, 10)
head(get_dimred(model), 5)
plot_dimred(model, label_milestones = T)

#### According to trajectory, Milestone number 4 is at the beginning of pseudotime
#### Milestone 4 will be selected as root for proper directionality
model <- model %>% add_root(root_milestone_id = "4")

#### create custom palette for plotting - our data contains 24 groups
c24 <- c(
  "#DC143C", #acinar
  "#FA8072", #ascl3
  "#A0522D",
  "#BA55D3",
  "#006400",
  "magenta4", #ID
  "goldenrod",
  "#0000CD", #ADULT SMGC
  "#6B8E23", # SD
  "#E6E6FA",  #E12
  "#66CDAA30",  #E12
  "#B0C4DE",  #E14
  "#66CDAA80",  #E14
  "#F4A460", #E16
  "#87CEEB",
  "#90EE90",
  "#F0E68C",
  "#F4A460", #P1
  "#DDA0DD",
  "#00FFFF", 
  "gold2",
  "#1E90FF")
pie(rep(1, 24), col = c24)
  
c24.2 <- c(
  "#E6E6FA",  #E12
  "#66CDAA30",  #E12
  "#B0C4DE",  #E14
  "#66CDAA80",  #E14
  "#F4A460", #E16
  "#87CEEB",
  "#90EE90",
  "#F0E68C",
  "#F4A460", #P1
  "#DDA0DD",
  "#8B0000", #doublepos
  "#00FFFF",
  "#32CD32",
  "gold2",
  "#1E90FF",
  "#DC143C", #acinar
  "#FA8072", #ascl3
  "#A0522D",
  "#BA55D3",
  "#006400",
  "magenta4", #ID
  "goldenrod",
  "#0000CD", #ADULT SMGC
  "#6B8E23")
pie(rep(1, 24), col = c24.2)

c15 <- c(
  "navy",
  "plum1",
  "purple",
  "chocolate2",
  "lightskyblue1",
  "goldenrod4",
  "orchid",
  "black",
  "greenyellow",
  "palegreen1",
  "royalblue3",
  "gold2",
  "firebrick1",
  "darkred",
  "magenta4"
  )
pie(rep(1, 15), col = c15)


#### PLOT PSEUDOTIME ###
plot_dimred(model, "pseudotime", pseudotime = calculate_pseudotime(model), grouping = stageids, color_density = "grouping", color_trajectory = "nearest", color_cells = "pseudotime" )
plot_graph(trajectory = model,arrow_length = unit(20, "pt"),label_milestones = T, grouping = dataset$grouping, color_cells = "grouping") + scale_y_reverse() + scale_color_manual(values = c24, aesthetics = "colour") + scale_fill_manual(values = c24, aesthetics = "fill")
plot_graph(trajectory = model, pseudotime = calculate_pseudotime(model), color_cells = "pseudotime", label_milestones = T, transition_size = 2,arrow_length = unit(20, "pt"))+ scale_y_reverse()
### plot trajectory
plot_dimred(model, grouping = dataset$grouping, color_density = "grouping", color_trajectory = "nearest", color_cells = "grouping", label_milestones = F, size_cells = 0.4, density_cutoff_label = 0, size_trajectory = 1,bw = 0.15, border_radius_percentage = 0.5) + scale_color_manual(values = c24, aesthetics = "colour") + scale_fill_manual(values = c24, aesthetics = "fill")
plot_dimred(model, grouping = groupsids, color_density = "grouping", color_trajectory = "nearest", color_cells = "grouping", label_milestones = F, size_cells = 0.4, density_cutoff_label = 0, size_trajectory = 1,bw = 0.15, border_radius_percentage = 0.5) + scale_color_manual(values = c15, aesthetics = "colour") + scale_fill_manual(values = c15, aesthetics = "fill")
plot_dimred(model, grouping = stageids, color_density = "grouping", color_trajectory = "nearest", color_cells = "grouping", label_milestones = F, size_cells = 0.4, density_cutoff_label = 0, size_trajectory = 1,bw = 0.15, border_radius_percentage = 0.5) 
plot_dimred(model, label_milestones = F, pseudotime = calculate_pseudotime(model), color_cells = "pseudotime", size_cells = 0.5, density_cutoff_label = 0, size_trajectory = 0.5,bw = 0.15, color_density = "grouping", grouping = dataset$grouping) 

### plot trajectory tree
plot_dendro(trajectory = model, size_cells = 2, color_cells = "pseudotime",pseudotime = calculate_pseudotime(model), y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link()
plot_dendro(trajectory = model, size_cells = 2, alpha_cells = 0.8, diag_offset = 0.005, grouping = dataset$grouping, color_cells = "grouping", y_offset = 0.6, color_milestones = "grouping", milestones = model$milestone_ids,border_radius_percentage = 0.1) + geom_edge_link() + scale_color_manual(values = c24, aesthetics = "colour") + scale_fill_manual(values = c24, aesthetics = "fill")
plot_dendro(trajectory = model, size_cells = 2, alpha_cells = 0.8, grouping = groupsids, color_cells = "grouping", y_offset = 0.6, color_milestones = "grouping", milestones = model$milestone_ids,border_radius_percentage = 0.1) + geom_edge_link() 
plot_dendro(trajectory = model, size_cells = 2, alpha_cells = 0.8, grouping = stageids, color_cells = "grouping", y_offset = 0.6, color_milestones = "grouping", milestones = model$milestone_ids,border_radius_percentage = 0.1) + geom_edge_link()

p<-plot_dendro(trajectory = model, size_cells = 2, alpha_cells = 0.8, grouping = stageids, color_cells = "grouping", y_offset = 0.6, color_milestones = "grouping", milestones = model$milestone_ids,border_radius_percentage = 0.1) + geom_edge_link()
g <- ggplot_build(p)
unique(g$data[[5]]["colour"])

c.stage <- c("#6495ED", "#20B2AA", "#DA70D6","#FA8072","#FFE4B5")
pie(rep(1, 5), col = c.stage)

plot_dendro(trajectory = model, size_cells = 2, alpha_cells = 0.8, grouping = stageids, color_cells = "grouping", y_offset = 0.6, color_milestones = "grouping", milestones = model$milestone_ids,border_radius_percentage = 0.1) + geom_edge_link() + scale_color_manual(values = c.stage, aesthetics = "colour") + scale_fill_manual(values = c.stage, aesthetics = "fill")
DimPlot(object = epithelium, reduction = "tsne", group.by = "stage", label = F, label.size = 6, pt.size = 1, cols = c.stage) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) 


### Final figures with color codes for consistency
## tree with stage and cell ID combined
pdf(file = "DYNO Trajectory Tree.pdf",width = 13,height = 9)
plot_dendro(trajectory = model, size_cells = 1, alpha_cells = 0.8, diag_offset = 0.005, grouping = dataset$grouping, color_cells = "grouping", y_offset = 0.6, color_milestones = "grouping", milestones = model$milestone_ids,border_radius_percentage = 0.1) + geom_edge_link() + scale_color_manual(values = c24, aesthetics = "colour") + scale_fill_manual(values = c24, aesthetics = "fill")
### with simple cell ID and stage
plot_dendro(trajectory = model, size_cells = 2, alpha_cells = 0.8, grouping = groupsids, color_cells = "grouping", y_offset = 0.6, color_milestones = "grouping", milestones = model$milestone_ids,border_radius_percentage = 0.1) + geom_edge_link() + scale_color_manual(values = c15, aesthetics = "colour") + scale_fill_manual(values = c15, aesthetics = "fill")
plot_dendro(trajectory = model, size_cells = 2, alpha_cells = 0.8, grouping = stageids, color_cells = "grouping", y_offset = 0.6, color_milestones = "grouping", milestones = model$milestone_ids,border_radius_percentage = 0.1) + geom_edge_link() + scale_color_manual(values = c("#DA70D6", "#FA8072","#FFE4B5", "#20B2AA",  "#6495ED"), aesthetics = "colour") + scale_fill_manual(values = c("#DA70D6", "#FA8072","#FFE4B5", "#20B2AA",  "#6495ED"), aesthetics = "fill")
dev.off()
### trajectory plots with complete ID and simple ID
pdf(file = "DYNO Trajectory.pdf",width = 8,height = 8)
plot_dimred(model, grouping = dataset$grouping, color_density = "none",alpha_cells = 0.3, color_trajectory = "nearest", color_cells = "grouping", label_milestones = F, size_cells = 2, density_cutoff_label = 0, size_trajectory = 1,bw = 0.15, border_radius_percentage = 0.5) + scale_color_manual(values = c24, aesthetics = "colour") + scale_fill_manual(values = c24, aesthetics = "fill")
plot_dimred(model, grouping = groupsids, color_density = "none", color_trajectory = "nearest", color_cells = "grouping", label_milestones = F, size_cells = 2, density_cutoff_label = 0, size_trajectory = 1,bw = 0.15, border_radius_percentage = 0.5) + scale_color_manual(values = c15, aesthetics = "colour") + scale_fill_manual(values = c15, aesthetics = "fill")
dev.off()

pdf(file = "DYNO Trajectory 2D graph.pdf",width = 8,height = 8)
plot_graph(trajectory = model,arrow_length = unit(20, "pt"),label_milestones = T, grouping = dataset$grouping, color_cells = "grouping") + scale_y_reverse() + scale_color_manual(values = c24, aesthetics = "colour") + scale_fill_manual(values = c24, aesthetics = "fill")
plot_graph(trajectory = model,arrow_length = unit(20, "pt"),label_milestones = T, grouping = groupsids, color_cells = "grouping") + scale_y_reverse() + scale_color_manual(values = c15, aesthetics = "colour") + scale_fill_manual(values = c15, aesthetics = "fill")
dev.off()

pdf(file = "tSNE epithelium by stage.pdf", width = 5.5,height = 4)
DimPlot(epithelium,pt.size = 1,group.by = "stage",cols = c("#6495ED", "#20B2AA", "#DA70D6","#FA8072","#FFE4B5")) + theme(axis.text = element_text(size = 16),axis.title = element_text(size = 16), legend.key.size = unit(x = 18,units = "pt"), legend.text = element_text(size = 16))
dev.off()

pdf(file = "tSNE epithelium by cell type.pdf", width = 8.5,height = 5)
DimPlot(epithelium,pt.size = 1,group.by = "Epi.CellType",cols = c15) + theme(axis.text = element_text(size = 16),axis.title = element_text(size = 16), legend.key.size = unit(x = 18,units = "pt"), legend.text = element_text(size = 16))
dev.off()

pdf(file = "tSNE epithelium by stage and cell type.pdf", width = 12,height = 5)
DimPlot(epithelium,pt.size = 1,group.by = "Epi.CellType.Stage",cols = c24) + theme(axis.text = element_text(size = 16),axis.title = element_text(size = 16), legend.key.size = unit(x = 18,units = "pt"), legend.text = element_text(size = 16))
dev.off()

pdf(file = "Dendrogram with expression of Transcription Factors.pdf", width = 13,height = 9)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Egr1", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Atf3", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Klf6", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Klf2", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Fos", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Jun", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Fosb", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Foxc1", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Foxi1", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Foxi2", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Bhlhe40", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Elf5", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Etv1", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Bhlha15", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Spdef", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Sox10", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Etv5", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Zbtb20", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Ehf", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Cebpd", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)

plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Cnn1", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Krt5", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)
plot_dendro(trajectory = model, size_cells = 2, feature_oi = "Gfra3", expression_source = dataset$expression, y_offset = 0.6, border_radius_percentage = 0.1) + geom_edge_link() +scale_fill_distiller(palette = "RdPu", aesthetics = c("colour", "fill"),direction = 1)


dev.off()
library(viridis)
########################
########################
########  tests  #######
########################
########################

### Calculate and plot overall chnages in features that change along the trajectory
overall_feature_importances <- dynfeature::calculate_overall_feature_importance(trajectory = model, expression_source = dataset$expression)

features <- overall_feature_importances %>% top_n(2000, importance) %>% pull(feature_id)
features.tfs <- features[features %in% Transcription.factors$Gene]
plot_onedim(add_root(model), grouping = dataset$grouping, color_cells = "grouping", size_cells = 2,alpha_cells = 0.8,label_milestones = T, orientation = -1)
plot_heatmap(model, expression_source = dataset, features_oi = features.tfs,grouping = groupsids,heatmap_type = "tiled", color_cells = "grouping") + scale_color_manual(values = c15, aesthetics = "colour") + scale_fill_manual(values = c15, aesthetics = "fill")

### Calculate and plot feature expression changes in each branch
branch_feature_importance <- calculate_branch_feature_importance(model, expression_source = dataset$expression,verbose = T, fi_method = )
###Features for branch 
branch_features.20 <- branch_feature_importance %>% filter(to == "20") %>% top_n(50, importance) %>% pull(feature_id)
plot_heatmap(model, expression_source = dataset, features_oi = branch_features.20,grouping = dataset$grouping,heatmap_type = "tiled")

### Calculate and plot feature expression changes AT bifurcation points
branch_feature_importance <- calculate_branch_feature_importance(model, expression_source = dataset$expression,verbose = T)
###Features for branch 
branching_milestone <- 6
branch6_feature_importance <- calculate_branching_point_feature_importance(model, milestones_oi = branching_milestone, expression_source = dataset$expression)
branch6.features <- branch_feature_importance %>% top_n(50, importance) %>% pull(feature_id)
plot_heatmap(model, expression_source = dataset, features_oi = unique(branch6.features),grouping = dataset$grouping,heatmap_type = "tiled")

