library(Seurat)
library(Hmisc)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
library(Matrix)

library(gridExtra)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

##### load file with merged epithelium 

epithelium <- readRDS(file = "~/2 - Integrate epithelial clusters/Merged epithelium from all stages (Seurat v3).rds")
epithelium <- SetIdent(epithelium, value = "Epi.CellType.Stage")
DimPlot(epithelium)

### load color scheme
c24 <- c("navy","plum1","purple", "chocolate2", "greenyellow", "palegreen1",
         "gold2","firebrick1","magenta4","gray","khaki1","black","yellow",
         "tan1","steelblue2","darkgreen","seagreen3","turquoise","lightskyblue1",
         "goldenrod4","orchid","royalblue3","darkgoldenrod3",
         "darkred")
pie(rep(1, 24), col = c24)

c24.2 <- c(
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
  "#8B0000", #doublepos
  "#00FFFF",
  "#32CD32",
  "gold2",
  "#1E90FF")
pie(rep(1, 24), col = c24.2)

c15 <- c( 
  "navy","plum1","purple","chocolate2",
  "lightskyblue1","goldenrod4","orchid","black","greenyellow","darkgreen",
  "royalblue3","gold2","firebrick1","darkred","magenta4"
)
pie(rep(1, 15), col = c15)

c9 <- c( 
  "khaki",
  "goldenrod3",
  "lightskyblue",
  "plum1",
  "royalblue3",
  "seagreen3",
  "darkred",
  "magenta4",
  "navy"
  
  )
pie(rep(1, 9), col = c9)

c12 <- c( 
  "khaki",
  "goldenrod3",
  "seagreen3",
  "lightskyblue1",
  "royalblue3",
  "darkgreen",
  "navy",
  "orchid",
  "darkred",
  "chocolate2",
  "black",
  "purple"
)
pie(rep(1, 12), col = c12)

#check that all metadata is correct and create plots
DimPlot(object = epithelium, reduction = "tsne", group.by = "seurat_clusters", label = T, pt.size = 1, label.size = 8, repel = F) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) + NoLegend()
DimPlot(object = epithelium, reduction = "tsne", group.by = "Epi.CellType", label = F, pt.size = 1, label.size = 5, repel = F, cols = c("navy", "hotpink", "darkred", "slateblue1", "slateblue3", "slateblue4", "firebrick2", "skyblue2", "goldenrod2", "purple", "firebrick", "forestgreen", "mediumpurple2", "pink2", "turquoise", "gray", "gray")) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) 
DimPlot(object = epithelium, reduction = "tsne", group.by = "Epi.CellType.Stage", label = F, pt.size = 1, label.size = 5, repel = F, cols = c("navy", "hotpink", "darkred","slateblue1", "goldenrod2", "purple2", "forestgreen", "pink2", "turquoise2", "lightskyblue1", "orange2", "lightskyblue3", "palevioletred3", "firebrick", "lightskyblue4", "orange", "green", "orangered3", "slateblue2", "darkorange2", "springgreen4", "orchid", "gray", "gray")) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) 
DimPlot(object = epithelium, reduction = "tsne", group.by = "stage", label = F, label.size = 6, pt.size = 1) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) 
DimPlot(object = epithelium, reduction = "tsne", group.by = "Epi.CellType", label = F, pt.size = 1, label.size = 5, repel = F, cols = c15) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) 
DimPlot(object = epithelium, reduction = "tsne", group.by = "Epi.CellType.Stage", label = F, pt.size = 1, label.size = 5, repel = F, cols = c24) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) 
DimPlot(object = epithelium, reduction = "tsne", group.by = "stage", label = F, label.size = 6, pt.size = 1, cols = ) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) 

epithelium <- SetIdent(epithelium, value = "Epi.CellType.Stage")
epithelium[["transition.labels"]] <- Idents(epithelium)
epithelium <- SetIdent(epithelium, value = "transition.labels")

DimPlot(object = epithelium, reduction = "tsne", group.by = "transition.labels", label = F, pt.size = 1, label.size = 5, repel = F) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) 

#### Load file from transcription factors
Transcription.factors <- read.delim("~/4 - Analysis of development transitions/Transcription factors.csv")

#########################################################
#########################################################
############  Analysis of duct populations ##############
epithelium <- SetIdent(epithelium, value = "Epi.CellType")
DimPlot(epithelium)
duct.epithelium <- subset(epithelium, idents = c("Krt19+ duct", "Basal duct", "Striated duct", "GCT", "Intercalated duct", "Ascl3", "Smgc+"))

plot1<-DimPlot(object = duct.epithelium, reduction = "tsne", group.by = "transition.labels", label = F, pt.size = 1, label.size = 5, repel = F,
        cols = c12) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18))

legend <- g_legend(plot1)

png(file = "tSNE duct populations.png", width = 1200,height = 500)
DimPlot(object = duct.epithelium, reduction = "tsne", group.by = "transition.labels", label = F, pt.size = 1,
        cols = c12) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) +NoLegend()
grid.arrange(legend,ncol=1, nrow=1)
dev.off()


epithelium <- SetIdent(epithelium, value = "transition.labels")

#########################################################
######### transitions from E12 to E14 ###################

ducte12.vs.ducte14 <- FindMarkers(object = epithelium, ident.1 = "E12_Krt19+ duct", ident.2 = "E14_Krt19+ duct", only.pos = F)
write.csv(ducte12.vs.ducte14, file = "E12 to E14 - duct Transition.csv", row.names = T)
ducte12.vs.ducte14.tfs <- ducte12.vs.ducte14[rownames(ducte12.vs.ducte14) %in% Transcription.factors$Gene, ]
write.csv(ducte12.vs.ducte14.tfs, file = "Pres Duct E12_vs_E14 - Transcription Factors.csv", row.names = T)
ducte12.vs.ducte14$gene <- rownames(ducte12.vs.ducte14)
ducte12.vs.ducte14.tfs$gene <- rownames(ducte12.vs.ducte14.tfs)

#### Generate violin plots
toptfs <- rbind((ducte12.vs.ducte14.tfs %>% top_n(5,avg_logFC)), (ducte12.vs.ducte14.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("E12_Krt19+ duct", "E14_Krt19+ duct"), cols = c13[1:2], features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend()
}
png("E12 to E14 Krt19 duct.png", width = 1500, height = 600,res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 

#########################################################
######### transitions from E14 to E16 ###################
#### Krt19 duct
ducte14.vs.k19e16 <- FindMarkers(object = epithelium, ident.1 = "E14_Krt19+ duct", ident.2 = "E16_Krt19+ duct", only.pos = F)
write.csv(ducte14.vs.k19e16, file = "E14 to E16 Krt19 duct transition.csv", row.names = T)
ducte14.vs.k19e16.tfs <- ducte14.vs.k19e16[rownames(ducte14.vs.k19e16) %in% Transcription.factors$Gene, ]
write.csv(ducte14.vs.k19e16.tfs, file = "E14_Krt19 to E16_Krt19 duct - Transcription Factors.csv", row.names = T)
ducte14.vs.k19e16$gene <- rownames(ducte14.vs.k19e16)
ducte14.vs.k19e16.tfs$gene <- rownames(ducte14.vs.k19e16.tfs)

#### Generate violin plots
toptfs <- rbind((ducte14.vs.k19e16.tfs %>% top_n(5,avg_logFC)), (ducte14.vs.k19e16.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("E14_Krt19+ duct", "E16_Krt19+ duct"), cols = c(c13[2],c13[4]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend()
}
png("E14 to E16 Krt19 duct.png", width = 1500, height = 600, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 

###### To E16 Basal duct
ducte14.vs.basalde16 <- FindMarkers(object = epithelium, ident.1 = "E14_Krt19+ duct", ident.2 = "E16_Basal duct", only.pos = F)
write.csv(ducte14.vs.basalde16, file = "E14 duct to E16 Basal duct.csv", row.names = T)
ducte14.vs.basalde16.tfs <- ducte14.vs.basalde16[rownames(ducte14.vs.basalde16) %in% Transcription.factors$Gene, ]
write.csv(ducte14.vs.basalde16.tfs, file = "E14 Krt19 duct to E16 Basal duct - Transcription Factors.csv", row.names = T)
ducte14.vs.basalde16$gene <- rownames(ducte14.vs.basalde16)
ducte14.vs.basalde16.tfs$gene <- rownames(ducte14.vs.basalde16.tfs)

#### Generate violin plots
toptfs <- rbind((ducte14.vs.basalde16.tfs %>% top_n(5,avg_logFC)), (ducte14.vs.basalde16.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("E14_Krt19+ duct", "E16_Basal duct"), cols = c(c13[2],c13[3]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("E14 Krt19 to E16 Basal duct.png", width = 1500, height = 600, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 


#########################################################
######### transitions from E16 to P1 ###################
#### Krt19 duct

k19ducte16.vs.k19ductp1 <- FindMarkers(object = epithelium, ident.1 = "E16_Krt19+ duct", ident.2 = c("P1_Krt19+ duct", "P1_Duct cluster 4"), only.pos = F)
write.csv(k19ducte16.vs.k19ductp1, file = "E16 to P1 Krt19 duct.csv", row.names = T)
k19ducte16.vs.k19ductp1.tfs <- k19ducte16.vs.k19ductp1[rownames(k19ducte16.vs.k19ductp1) %in% Transcription.factors$Gene, ]
write.csv(k19ducte16.vs.k19ductp1.tfs, file = "E16 Krt19 duct to P1 Krt19 duct - Transcription Factors.csv", row.names = T)
k19ducte16.vs.k19ductp1$gene <- rownames(k19ducte16.vs.k19ductp1)
k19ducte16.vs.k19ductp1.tfs$gene <- rownames(k19ducte16.vs.k19ductp1.tfs)

#### Generate violin plots
toptfs <- rbind((k19ducte16.vs.k19ductp1.tfs %>% top_n(5,avg_logFC)), (k19ducte16.vs.k19ductp1.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("E16_Krt19+ duct", "P1_Krt19+ duct", "P1_Duct cluster 4"), cols = c(c13[4:5], c13[7]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("E16 Krt19 to P1 Krt19 duct.png", width = 1500, height = 600, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 

### Basal duct
basalde16.vs.basaldp1 <- FindMarkers(object = epithelium, ident.1 = "E16_Basal duct", ident.2 = "P1_Basal duct", only.pos = F)
write.csv(basalde16.vs.basaldp1, file = "E16 vs P1 Basal duct.csv", row.names = T)
basalde16.vs.basaldp1.tfs <- basalde16.vs.basaldp1[rownames(basalde16.vs.basaldp1) %in% Transcription.factors$Gene, ]
write.csv(basalde16.vs.basaldp1.tfs, file = "E16_Basal duct to P1 Basal duct - Transcription Factors.csv", row.names = T)
basalde16.vs.basaldp1$gene <- rownames(basalde16.vs.basaldp1)
basalde16.vs.basaldp1.tfs$gene <- rownames(basalde16.vs.basaldp1.tfs)

toptfs <- rbind((basalde16.vs.basaldp1.tfs %>% top_n(5,avg_logFC)), (basalde16.vs.basaldp1.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("E16_Basal duct", "P1_Basal duct"), cols = c(c13[3],c13[6]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("E16 Basal duct to P1 Basal duct.png", width = 1500, height = 600, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 

###Postnatal duct differentiation (P1 - to Adult) ####
### Basal duct
basaldp1.vs.basaldAdult <- FindMarkers(object = epithelium, ident.1 = "P1_Basal duct", ident.2 = "Adult_Basal duct", only.pos = F)
write.csv(basaldp1.vs.basaldAdult, file = "P1 to Adult Basal duct.csv", row.names = T)
basaldp1.vs.basaldAdult.tfs <- basaldp1.vs.basaldAdult[rownames(basaldp1.vs.basaldAdult) %in% Transcription.factors$Gene, ]
write.csv(basaldp1.vs.basaldAdult.tfs, file = "P1_Basal duct to Adult Basal duct - Transcription Factors.csv", row.names = T)
basaldp1.vs.basaldAdult$gene <- rownames(basaldp1.vs.basaldAdult)
basaldp1.vs.basaldAdult.tfs$gene <- rownames(basaldp1.vs.basaldAdult.tfs)

toptfs <- rbind((basaldp1.vs.basaldAdult.tfs %>% top_n(5,avg_logFC)), (basaldp1.vs.basaldAdult.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("P1_Basal duct", "Adult_Basal duct"), cols = c(c13[6],c13[12]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("P1 Basal duct to Adult Basal duct.png", width = 1500, height = 600, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 

#### Intercalated duct
k19ductp1.vs.idp30 <- FindMarkers(object = epithelium, ident.1 = c("P1_Krt19+ duct"), ident.2 = "Adult_Intercalated duct", only.pos = F)
write.csv(k19ductp1.vs.idp30, file = "P1 duct to Intercalated duct.csv", row.names = T)
k19ductp1.vs.idp30.tfs <- k19ductp1.vs.idp30[rownames(k19ductp1.vs.idp30) %in% Transcription.factors$Gene, ]
write.csv(k19ductp1.vs.idp30.tfs, file = "P1 duct to Intercalated duct - Transcription Factors.csv", row.names = T)
k19ductp1.vs.idp30$gene <- rownames(k19ductp1.vs.idp30)
k19ductp1.vs.idp30.tfs$gene <- rownames(k19ductp1.vs.idp30.tfs)

toptfs <- rbind((k19ductp1.vs.idp30.tfs %>% top_n(5,avg_logFC)), (k19ductp1.vs.idp30.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("P1_Krt19+ duct", "Adult_Intercalated duct" ), cols = c(c13[5],c13[10]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) + NoLegend()
}
png("P1 Duct to adult Intercalated duct.png", width = 1500, height = 600, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 

#### Ascl3+ cells
p1duct.toascl3 <- FindMarkers(object = epithelium, ident.1 = c("P1_Krt19+ duct"), ident.2 = "Adult_Ascl3", only.pos = F)
write.csv(p1duct.toascl3, file = "P1 duct to adult Ascl3.csv", row.names = T)
p1duct.toascl3.tfs <- p1duct.toascl3[rownames(p1duct.toascl3) %in% Transcription.factors$Gene, ]
write.csv(p1duct.toascl3.tfs, file = "P1 duct to adult Ascl3 - Transcription Factors.csv", row.names = T)
p1duct.toascl3$gene <- rownames(p1duct.toascl3)
p1duct.toascl3.tfs$gene <- rownames(p1duct.toascl3.tfs)

toptfs <- rbind((p1duct.toascl3.tfs %>% top_n(5,avg_logFC)), (p1duct.toascl3.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("P1_Krt19+ duct", "Adult_Ascl3" ), cols = c(c13[5],c13[9]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("P1 Duct to adult Ascl3.png", width = 1500, height = 600, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 

#### Striated duct
sdp1.to.sdp30 <- FindMarkers(object = epithelium, ident.1 = c("P1_Krt19+ duct"), ident.2 = "Adult_Striated duct", only.pos = F)
sdp1.to.sdp30 <- sdp1.to.sdp30[sdp1.to.sdp30$p_val_adj <0.05, ]
write.csv(sdp1.to.sdp30, file = "P1 duct to adult Striated duct.csv", row.names = T)
sdp1.to.sdp30.tfs <- sdp1.to.sdp30[rownames(sdp1.to.sdp30) %in% Transcription.factors$Gene, ]
write.csv(sdp1.to.sdp30.tfs, file = "P1 duct to adult Striated duct - Transcription Factors.csv", row.names = T)
sdp1.to.sdp30$gene <- rownames(sdp1.to.sdp30)
sdp1.to.sdp30.tfs$gene <- rownames(sdp1.to.sdp30.tfs)

toptfs <- sdp1.to.sdp30.tfs %>% top_n(5,avg_logFC)
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("P1_Krt19+ duct",  "Adult_Striated duct" ), cols = c(c13[5],c13[13]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("P1 Duct to adult Striated duct.png", width = 1500, height = 350, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 

#### gct

sdp1.to.gctp30 <- FindMarkers(object = epithelium, ident.1 = c("P1_Krt19+ duct"), ident.2 = "Adult_GCT", only.pos = F)
write.csv(sdp1.to.gctp30, file = "P1 duct to GCT.csv", row.names = T)
sdp1.to.gctp30.tfs <- sdp1.to.gctp30[rownames(sdp1.to.gctp30) %in% Transcription.factors$Gene, ]
write.csv(sdp1.to.gctp30.tfs, file = "P1 duct to GCT - Transcription Factors.csv", row.names = T)
sdp1.to.gctp30$gene <- rownames(sdp1.to.gctp30)
sdp1.to.gctp30.tfs$gene <- rownames(sdp1.to.gctp30.tfs)

toptfs <- rbind((sdp1.to.gctp30.tfs %>% top_n(5,avg_logFC)), (sdp1.to.gctp30.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("P1_Krt19+ duct", "Adult_GCT" ), cols = c(c13[5],c13[11]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("P1 Duct to adult GCT.png", width = 1500, height = 600, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 


#### Smgc+ cells
krt19.vs.smgcp30 <- FindMarkers(object = epithelium, ident.1 = c("P1_Krt19+ duct"), ident.2 = "Adult_Smgc+", only.pos = F)
write.csv(krt19.vs.smgcp30, file = "P1 duct to Adult Smgc.csv", row.names = T)
krt19.vs.smgcp30.tfs <- krt19.vs.smgcp30[rownames(krt19.vs.smgcp30) %in% Transcription.factors$Gene, ]
write.csv(krt19.vs.smgcp30.tfs, file = "P1 duct to Adult Smgc - Transcription Factors.csv", row.names = T)
krt19.vs.smgcp30$gene <- rownames(krt19.vs.smgcp30)
krt19.vs.smgcp30.tfs$gene <- rownames(krt19.vs.smgcp30.tfs)

toptfs <- rbind((krt19.vs.smgcp30.tfs %>% top_n(5,avg_logFC)), (krt19.vs.smgcp30.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("P1_Krt19+ duct", "Adult_Smgc+" ), cols = c(c13[5],c13[8]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("P1 Duct to adult Smgc.png", width = 1500, height = 600, res=300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 


Smgc.vs.smgcp30 <- FindMarkers(object = epithelium, ident.1 = c("P1_Smgc+ Proacinar"), ident.2 = "Adult_Smgc+", only.pos = F)
write.csv(Smgc.vs.smgcp30, file = "P1 Smgc to Adult Smgc.csv", row.names = T)
Smgc.vs.smgcp30.tfs <- Smgc.vs.smgcp30[rownames(Smgc.vs.smgcp30) %in% Transcription.factors$Gene, ]
write.csv(Smgc.vs.smgcp30.tfs, file = "P1 smgc to Adult Smgc - Transcription Factors.csv", row.names = T)
Smgc.vs.smgcp30$gene <- rownames(Smgc.vs.smgcp30)
Smgc.vs.smgcp30.tfs$gene <- rownames(Smgc.vs.smgcp30.tfs)

toptfs <- rbind((Smgc.vs.smgcp30.tfs %>% top_n(5,avg_logFC)), (Smgc.vs.smgcp30.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("P1_Smgc+ Proacinar", "Adult_Smgc+" ), cols = c("pink2",c13[8]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend()
}
png("P1 smgc to adult Smgc.png", width = 1500, height = 600, res=300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 


#### Basal duct to SD and GCT

epithelium <- SetIdent(epithelium, value = "transition.labels")
DimPlot(epithelium)
DimPlot(object = duct.epithelium, reduction = "tsne", group.by = "transition.labels", label = F, pt.size = 1, label.size = 5, repel = F,
               cols = c13) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18))

bd.to.sd <- FindMarkers(object = epithelium, ident.1 = "Adult_Basal duct", ident.2 = "Adult_Striated duct", only.pos = F)
bd.to.sd <- bd.to.sd[bd.to.sd$p_val_adj <0.05, ]
write.csv(bd.to.sd, file = "Basal duct to striated duct.csv", row.names = T)
bd.to.sd.tfs <- bd.to.sd[rownames(bd.to.sd) %in% Transcription.factors$Gene, ]
write.csv(bd.to.sd.tfs, file = "Basal duct to striated duct - Transcription Factors.csv", row.names = T)
bd.to.sd$gene <- rownames(bd.to.sd)
bd.to.sd.tfs$gene <- rownames(bd.to.sd.tfs)

epithelium <- SetIdent(epithelium, value = "Epi.CellType")
toptfs <- rbind((bd.to.sd.tfs %>% top_n(5,avg_logFC)), (bd.to.sd.tfs %>% top_n(-1,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = T),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("Basal duct", "Striated duct"), cols = c(c13[12],c13[13]), features = toptfs$gene[i], pt.size = 0, group.by = "Epi.CellType") + theme(axis.title = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("Basal Duct to Striated duct.png", width = 1500, height = 600, res=300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 

epithelium <- SetIdent(epithelium, value = "transition.labels")
bd.to.gct <- FindMarkers(object = epithelium, ident.1 = "Adult_Basal duct", ident.2 = "Adult_GCT", only.pos = F)
bd.to.gct <- bd.to.gct[bd.to.gct$p_val_adj <0.05, ]
write.csv(bd.to.gct, file = "Basal duct to GCT.csv", row.names = T)
bd.to.gct.tfs <- bd.to.gct[rownames(bd.to.gct) %in% Transcription.factors$Gene, ]
write.csv(bd.to.gct.tfs, file = "Basal duct to GCT - Transcription Factors.csv", row.names = T)
bd.to.gct$gene <- rownames(bd.to.gct)
bd.to.gct.tfs$gene <- rownames(bd.to.gct.tfs)

epithelium <- SetIdent(epithelium, value = "Epi.CellType")
toptfs <- rbind((bd.to.gct.tfs %>% top_n(5,avg_logFC)), (bd.to.gct.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("Basal duct", "GCT"), cols = c(c13[12],c13[11]), features = toptfs$gene[i], pt.size = 0, group.by = "Epi.CellType") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("Basal Duct to GCT.png", width = 1500, height = 600, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 


##END BUD/ACINAR POPULATIONS##
##############################
#Create a temporary URD object just to create TSNE plots for visualization -- requires previously generated URD object#
epithelium <- SetIdent(epithelium, value = "Epi.CellType")
DimPlot(epithelium)
acinar.epithelium <- subset(epithelium, idents = c("End bud", "Smgc+ Proacinar", "Bpifa2+ Proacinar", "Bpifa2+ Smgc+ Proacinar", "Acinar", "Bpifa2+ Acinar"))
plot1<-DimPlot(object = acinar.epithelium, reduction = "tsne", group.by = "transition.labels", label = F, pt.size = 1, label.size = 5, repel = F,
        cols = c9) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) 
legend <- g_legend(plot1)

png(file = "tSNE acinar populations.png", width = 1500.5,height = 5)
DimPlot(object = acinar.epithelium, reduction = "tsne", group.by = "transition.labels", label = F, pt.size = 1, label.size = 5, repel = F,
        cols = c9) + theme(axis.text = element_text(size=18), axis.title = element_text(size = 18)) +NoLegend()
grid.arrange(legend,ncol=1, nrow=1)
dev.off()


epithelium <- SetIdent(epithelium, value = "transition.labels")
###E12 TO E14####

endbud12.vs.endbude14 <- FindMarkers(object = epithelium, ident.1 = "E12_End bud", ident.2 = "E14_End bud", only.pos = F)
write.csv(endbud12.vs.endbude14, file = "E12 End bud to E14.csv", row.names = T)
endbud12.vs.endbude14.tfs <- endbud12.vs.endbude14[rownames(endbud12.vs.endbude14) %in% Transcription.factors$Gene, ]
write.csv(endbud12.vs.endbude14.tfs, file = "E12 End bud to E14 - Transcription Factors.csv", row.names = T)
endbud12.vs.endbude14$gene <- rownames(endbud12.vs.endbude14)
endbud12.vs.endbude14.tfs$gene <- rownames(endbud12.vs.endbude14.tfs)

toptfs <- rbind((endbud12.vs.endbude14.tfs %>% top_n(5,avg_logFC)), (endbud12.vs.endbude14.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(acinar.epithelium, idents = c("E12_End bud", "E14_End bud"), cols = c9[1:2], features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend()
}
png("E12 End bud to E14.png", width = 1500, height = 600, res=300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 


###E14 TO E16####
endbude14.vs.endbude16 <- FindMarkers(object = epithelium, ident.1 = "E14_End bud", ident.2 = "E16_End bud", only.pos = F)
write.csv(endbude14.vs.endbude16, file = "E14 End bud  to E16.csv", row.names = T)
endbude14.vs.endbude16.tfs <- endbude14.vs.endbude16[rownames(endbude14.vs.endbude16) %in% Transcription.factors$Gene, ]
write.csv(endbude14.vs.endbude16.tfs, file = "E14 End bud  to E16 - Transcription Factors.csv", row.names = T)
endbude14.vs.endbude16$gene <- rownames(endbude14.vs.endbude16)
endbude14.vs.endbude16.tfs$gene <- rownames(endbude14.vs.endbude16.tfs)

toptfs <- rbind((endbude14.vs.endbude16.tfs %>% top_n(5,avg_logFC)), (endbude14.vs.endbude16.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("E14_End bud", "E16_End bud"), cols = c9[2:3], features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend()
}
png("E14 End bud to E16.png", width = 1500, height = 600, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 


###E16 TO P1####
endbude16.vs.bpifa2 <- FindMarkers(object = epithelium, ident.1 = "E16_End bud", ident.2 = c("P1_Bpifa2+ Proacinar"), only.pos = F)
write.csv(endbude16.vs.bpifa2, file = "E16 End bud to P1 Bpifa2 Proacinar.csv", row.names = T)
endbude16.vs.bpifa2.tfs <- endbude16.vs.bpifa2[rownames(endbude16.vs.bpifa2) %in% Transcription.factors$Gene, ]
write.csv(endbude16.vs.bpifa2.tfs, file = "E16 End bud to P1 Bpifa2 Proacinar - Transcription Factors.csv", row.names = T)
endbude16.vs.bpifa2$gene <- rownames(endbude16.vs.bpifa2)
endbude16.vs.bpifa2.tfs$gene <- rownames(endbude16.vs.bpifa2.tfs)

toptfs <- rbind((endbude16.vs.bpifa2.tfs %>% top_n(5,avg_logFC)), (endbude16.vs.bpifa2.tfs %>% top_n(-4,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("E16_End bud", "P1_Bpifa2+ Proacinar", "P1_Bpifa2+ Smgc+ Proacinar"), cols = c(c9[3],c9[5:6]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("E16 End bud to P1 Bpifa2.png", width = 1500, height = 600, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 


endbude16.vs.smgc <- FindMarkers(object = epithelium, ident.1 = "E16_End bud", ident.2 = c("P1_Smgc+ Proacinar"), only.pos = F)
write.csv(endbude16.vs.smgc, file = "E16 End bud to P1 smgc Proacinar.csv", row.names = T)
endbude16.vs.smgc.tfs <- endbude16.vs.smgc[rownames(endbude16.vs.smgc) %in% Transcription.factors$Gene, ]
write.csv(endbude16.vs.smgc.tfs, file = "E16 End bud to P1 smgc Proacinar - Transcription Factors.csv", row.names = T)
endbude16.vs.smgc$gene <- rownames(endbude16.vs.smgc)
endbude16.vs.smgc.tfs$gene <- rownames(endbude16.vs.smgc.tfs)

toptfs <- rbind((endbude16.vs.smgc.tfs %>% top_n(5,avg_logFC)), (endbude16.vs.smgc.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("E16_End bud", "P1_Smgc+ Proacinar", "P1_Bpifa2+ Smgc+ Proacinar"), cols = c(c9[3:4],c9[6]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("E16 End bud to P1 smgc.png", width = 1500, height = 600, res = 300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off() 


###P1 to P30####
p1.bpifa2.to.adult <- FindMarkers(object = epithelium, ident.1 = "P1_Bpifa2+ Proacinar", ident.2 = "Adult_Bpifa2+ Acinar", only.pos = F)
p1.bpifa2.to.adult <- p1.bpifa2.to.adult[p1.bpifa2.to.adult$p_val_adj<0.05, ]
write.csv(p1.bpifa2.to.adult, file = "P1 Bpifa2 Proacinar to Adult Bpifa2.csv", row.names = T)
p1.bpifa2.to.adult.tfs <- p1.bpifa2.to.adult[rownames(p1.bpifa2.to.adult) %in% Transcription.factors$Gene, ]
write.csv(p1.bpifa2.to.adult.tfs, file = "P1 Bpifa2 Proacinar to Adult Bpifa2 - Transcription Factors.csv", row.names = T)
p1.bpifa2.to.adult$gene <- rownames(p1.bpifa2.to.adult)
p1.bpifa2.to.adult.tfs$gene <- rownames(p1.bpifa2.to.adult.tfs)

toptfs <- rbind((p1.bpifa2.to.adult.tfs %>% top_n(5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("P1_Bpifa2+ Proacinar", "Adult_Bpifa2+ Acinar"), cols = c(c9[5],c9[8]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("P1 Bpifa2 to Adult Bpifa2.png", width = 1500, height = 350, res=300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off()

##smgc
p1.smgc.to.adult <- FindMarkers(object = epithelium, ident.1 = "P1_Smgc+ Proacinar", ident.2 = "Adult_Smgc+", only.pos = F)
p1.smgc.to.adult <- p1.smgc.to.adult[p1.smgc.to.adult$p_val_adj <0.05,]
write.csv(p1.smgc.to.adult, file = "P1 Smgc Proacinar to Adult Smgc.csv", row.names = T)
p1.smgc.to.adult.tfs <- p1.smgc.to.adult[rownames(p1.smgc.to.adult) %in% Transcription.factors$Gene, ]
write.csv(p1.smgc.to.adult.tfs, file = "P1 Smgc Proacinar to Adult Smgc - Transcription Factors.csv", row.names = T)
p1.smgc.to.adult$gene <- rownames(p1.smgc.to.adult)
p1.smgc.to.adult.tfs$gene <- rownames(p1.smgc.to.adult.tfs)

toptfs <- rbind((p1.smgc.to.adult.tfs %>% top_n(5,avg_logFC)), (p1.smgc.to.adult.tfs %>% top_n(-5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("P1_Smgc+ Proacinar", "Adult_Smgc+"), cols = c(c9[4],c9[9]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("P1 Smgc to Adult Smgc.png", width = 1500, height = 600, res=300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off()


p1.proacinar.to.acinar <- FindMarkers(object = epithelium, ident.1 = c("P1_Smgc+ Proacinar", "P1_Bpifa2+ Proacinar"), ident.2 = "Adult_Acinar", only.pos = F)
write.csv(p1.proacinar.to.acinar, file = "P1 Proacinar to Adult Acinar.csv", row.names = T)
p1.proacinar.to.acinar.tfs <- p1.proacinar.to.acinar[rownames(p1.proacinar.to.acinar) %in% Transcription.factors$Gene, ]
write.csv(p1.proacinar.to.acinar.tfs, file = "P1 Proacinar to Adult Acinar - Transcription Factors.csv", row.names = T)
p1.proacinar.to.acinar$gene <- rownames(p1.proacinar.to.acinar)
p1.proacinar.to.acinar.tfs$gene <- rownames(p1.proacinar.to.acinar.tfs)

toptfs <- rbind((p1.proacinar.to.acinar.tfs %>% top_n(5,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]

plot.list <- list()
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("P1_Smgc+ Proacinar", "P1_Bpifa2+ Proacinar", "Adult_Acinar"), cols = c(c9[4:5],c9[7]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("P1 proacinar to Adult Acinar.png", width = 1500, height = 350, res=300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off()

 
p1.k19duct.to.acinar <- FindMarkers(object = epithelium, ident.1 = c("P1_Krt19+ duct"), ident.2 = "Adult_Acinar", only.pos = F)
write.csv(p1.k19duct.to.acinar, file = "P1 K19 duct to Adult Acinar.csv", row.names = T)
p1.k19duct.to.acinar.tfs <- p1.k19duct.to.acinar[rownames(p1.k19duct.to.acinar) %in% Transcription.factors$Gene, ]
write.csv(p1.k19duct.to.acinar.tfs, file = "P1 K19 duct to Adult Acinar - Transcription Factors.csv", row.names = T)
p1.k19duct.to.acinar$gene <- rownames(p1.k19duct.to.acinar)
p1.k19duct.to.acinar.tfs$gene <- rownames(p1.k19duct.to.acinar.tfs)

toptfs <- rbind((p1.k19duct.to.acinar.tfs %>% top_n(5,avg_logFC)), (p1.k19duct.to.acinar.tfs %>% top_n(-2,avg_logFC)))
toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]

plot.list <- list()        
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("P1_Krt19+ duct", "Adult_Acinar"), cols = c("royalblue3", c9[7]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
}
png("P1 k19 duct to Adult Acinar.png", width = 1500, height = 600, res=300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off()


####################### EXPLORATORY ANALYSIS ############################
##### POTENTIAL TRANSITIONS IN ADULT POPULATIONS ######
DimPlot(epithelium)

#### Acinar to/from Bpifa2
acinar.vs.bpifa2 <- FindMarkers(object = epithelium, ident.1 = c("Adult_Acinar"), ident.2 = "Adult_Bpifa2+ Acinar", only.pos = F)
write.csv(acinar.vs.bpifa2, file = "Transdiff potential - Acinar-Bpifa2.csv", row.names = T)
acinar.vs.bpifa2.tfs <- acinar.vs.bpifa2[rownames(acinar.vs.bpifa2) %in% Transcription.factors$Gene, ]
write.csv(acinar.vs.bpifa2.tfs, file = "Transdiff potential - Acinar-Bpifa2 - Transcription Factors.csv", row.names = T)
acinar.vs.bpifa2$gene <- rownames(acinar.vs.bpifa2)
acinar.vs.bpifa2.tfs$gene <- rownames(acinar.vs.bpifa2.tfs)

(acinar.vs.bpifa2.tfs %>% top_n(-5,avg_logFC))

#toptfs <- rbind((acinar.vs.bpifa2.tfs %>% top_n(5,avg_logFC)), (acinar.vs.bpifa2.tfs %>% top_n(-5,avg_logFC)))
#toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]
#plot.list <- list()        
#for(i in 1:length(toptfs$gene)){
#  plot.list[[i]] <- VlnPlot(epithelium, idents = c("P1_Krt19+ duct", "Adult_Acinar"), cols = c("royalblue3", c9[7]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend() 
#}
#png("P1 k19 duct to Adult Acinar.png", width = 1500, height = 600, res=300) 
#CombinePlots(plots = plot.list, ncol = 5)
#dev.off()


#### Acinar to/from Smgc
acinar.vs.smgc <- FindMarkers(object = epithelium, ident.1 = c("Adult_Acinar"), ident.2 = "Adult_Smgc+", only.pos = F)
write.csv(acinar.vs.smgc, file = "Transdiff potential - Acinar-Smgc.csv", row.names = T)
acinar.vs.smgc.tfs <- acinar.vs.smgc[rownames(acinar.vs.smgc) %in% Transcription.factors$Gene, ]
write.csv(acinar.vs.smgc.tfs, file = "Transdiff potential - Acinar-smgc - Transcription Factors.csv", row.names = T)
acinar.vs.smgc$gene <- rownames(acinar.vs.smgc)
acinar.vs.smgc.tfs$gene <- rownames(acinar.vs.smgc.tfs)
toptfs <- rbind((acinar.vs.smgc.tfs %>% top_n(5,avg_logFC)), (acinar.vs.smgc.tfs %>% top_n(-5,avg_logFC)))

#### Acinar to/from ID
acinar.vs.ID <- FindMarkers(object = epithelium, ident.1 = c("Adult_Acinar"), ident.2 = "Adult_Intercalated duct", only.pos = F)
write.csv(acinar.vs.ID, file = "Transdiff potential - Acinar-Intercalated duct.csv", row.names = T)
acinar.vs.ID.tfs <- acinar.vs.ID[rownames(acinar.vs.ID) %in% Transcription.factors$Gene, ]
write.csv(acinar.vs.ID.tfs, file = "Transdiff potential - Acinar-ID - Transcription Factors.csv", row.names = T)
acinar.vs.ID$gene <- rownames(acinar.vs.ID)
acinar.vs.ID.tfs$gene <- rownames(acinar.vs.ID.tfs)
toptfs <- rbind((acinar.vs.ID.tfs %>% top_n(5,avg_logFC)), (acinar.vs.ID.tfs %>% top_n(-5,avg_logFC)))

#### Acinar to/from Basal duct
acinar.vs.BD <- FindMarkers(object = epithelium, ident.1 = c("Adult_Acinar"), ident.2 = "Adult_Basal duct", only.pos = F)
write.csv(acinar.vs.BD, file = "Transdiff potential - Acinar-Basal duct.csv", row.names = T)
acinar.vs.BD.tfs <- acinar.vs.BD[rownames(acinar.vs.BD) %in% Transcription.factors$Gene, ]
write.csv(acinar.vs.BD.tfs, file = "Transdiff potential - Acinar-BD - Transcription Factors.csv", row.names = T)
acinar.vs.BD$gene <- rownames(acinar.vs.BD)
acinar.vs.BD.tfs$gene <- rownames(acinar.vs.BD.tfs)
toptfs <- rbind((acinar.vs.BD.tfs %>% top_n(2,avg_logFC)), (acinar.vs.BD.tfs %>% top_n(-5,avg_logFC)))

toptfs <- toptfs[order(toptfs$avg_logFC, decreasing = F),]

plot.list <- list()        
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("Adult_Basal duct", "Adult_Acinar"), cols = c("black", c9[7]), features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels", sort = "decreasing") + theme(axis.title = element_blank(), axis.text.x =element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend()
}
png("Basal duct to Adult Acinar.png", width = 1500, height = 600, res=300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off()




#### Acinar to/from GCT
acinar.vs.GCT <- FindMarkers(object = epithelium, ident.1 = c("Adult_Acinar"), ident.2 = "Adult_GCT", only.pos = F)
write.csv(acinar.vs.GCT, file = "Transdiff potential - Acinar-GCT.csv", row.names = T)
acinar.vs.GCT.tfs <- acinar.vs.GCT[rownames(acinar.vs.GCT) %in% Transcription.factors$Gene, ]
write.csv(acinar.vs.GCT.tfs, file = "Transdiff potential - Acinar-GCT - Transcription Factors.csv", row.names = T)
acinar.vs.GCT$gene <- rownames(acinar.vs.GCT)
acinar.vs.GCT.tfs$gene <- rownames(acinar.vs.GCT.tfs)
toptfs <- rbind((acinar.vs.GCT.tfs %>% top_n(5,avg_logFC)), (acinar.vs.GCT.tfs %>% top_n(-5,avg_logFC)))


#### Acinar to/from SD
epithelium <- SetIdent(epithelium, value = "Epi.CellType.Stage")
acinar.vs.SD <- FindMarkers(object = epithelium, ident.1 = c("Adult_Acinar"), ident.2 = "Adult_Striated duct", only.pos = F, min.pct = 0.25)
acinar.vs.SD <- acinar.vs.SD[acinar.vs.SD$p_val_adj<0.05, ]
write.csv(acinar.vs.SD, file = "Transdiff potential - Acinar-SD.csv", row.names = T)
acinar.vs.SD.tfs <- acinar.vs.SD[rownames(acinar.vs.SD) %in% Transcription.factors$Gene, ]
write.csv(acinar.vs.SD.tfs, file = "Transdiff potential - Acinar-SD - Transcription Factors.csv", row.names = T)
acinar.vs.SD$gene <- rownames(acinar.vs.SD)
acinar.vs.SD.tfs$gene <- rownames(acinar.vs.SD.tfs)

toptfs <- acinar.vs.SD.tfs[c(1:7),]
toptfs <-  toptfs[order(toptfs$avg_logFC, decreasing = F),]

plot.list <- list()        
for(i in 1:length(toptfs$gene)){
  plot.list[[i]] <- VlnPlot(epithelium, idents = c("Adult_Striated duct", "Adult_Acinar"), cols = c("#00CC00", "#CC0000"),sort = "decreasing", features = toptfs$gene[i], pt.size = 0, group.by = "transition.labels") + theme(axis.title = element_blank(), axis.text.x  = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(), title = element_text(size = 10)) +NoLegend()
}
png("Adult striated duct to Adult Acinar.png", width = 1500, height = 600, res=300) 
CombinePlots(plots = plot.list, ncol = 5)
dev.off()





#### Acinar to/from Ascl3
acinar.vs.Ascl3 <- FindMarkers(object = epithelium, ident.1 = c("Adult_Acinar"), ident.2 = "Adult_Ascl3", only.pos = F)
write.csv(acinar.vs.Ascl3, file = "Transdiff potential - Acinar-Ascl3.csv", row.names = T)
acinar.vs.Ascl3.tfs <- acinar.vs.Ascl3[rownames(acinar.vs.Ascl3) %in% Transcription.factors$Gene, ]
write.csv(acinar.vs.Ascl3.tfs, file = "Transdiff potential - Acinar-Ascl3 - Transcription Factors.csv", row.names = T)
acinar.vs.Ascl3$gene <- rownames(acinar.vs.Ascl3)
acinar.vs.Ascl3.tfs$gene <- rownames(acinar.vs.Ascl3.tfs)
toptfs <- rbind((acinar.vs.Ascl3.tfs %>% top_n(5,avg_logFC)), (acinar.vs.Ascl3.tfs %>% top_n(-5,avg_logFC)))

#### Intercalated duct to/from SMGC 
ID.to.Smgc <- FindMarkers(object = epithelium, ident.1 = c("Adult_Intercalated duct"), ident.2 = "Adult_Smgc+", only.pos = F)
write.csv(ID.to.Smgc, file = "Transdiff potential - Intercalated duct to Smgc.csv", row.names = T)
ID.to.Smgc.tfs <- ID.to.Smgc[rownames(ID.to.Smgc) %in% Transcription.factors$Gene, ]
write.csv(ID.to.Smgc.tfs, file = "Transdiff potential - Intercalated duct to Smgc - Transcription Factors.csv", row.names = T)
ID.to.Smgc$gene <- rownames(ID.to.Smgc)
ID.to.Smgc.tfs$gene <- rownames(ID.to.Smgc.tfs)
toptfs <- rbind((ID.to.Smgc.tfs %>% top_n(5,avg_logFC)), (ID.to.Smgc.tfs %>% top_n(-5,avg_logFC)))

toptfs <-  toptfs[order(toptfs$avg_logFC, decreasing = F),]








#### Intercalated duct to/from Bpifa2 
ID.to.Bpifa2 <- FindMarkers(object = epithelium, ident.1 = c("Adult_Intercalated duct"), ident.2 = "Adult_Bpifa2+ Acinar", only.pos = F)
write.csv(ID.to.Bpifa2, file = "Transdiff potential - Intercalated duct to Bpifa2.csv", row.names = T)
ID.to.Bpifa2.tfs <- ID.to.Bpifa2[rownames(ID.to.Bpifa2) %in% Transcription.factors$Gene, ]
write.csv(ID.to.Bpifa2.tfs, file = "Transdiff potential - Intercalated duct to Bpifa2 - Transcription Factors.csv", row.names = T)
ID.to.Bpifa2$gene <- rownames(ID.to.Bpifa2)
ID.to.Bpifa2.tfs$gene <- rownames(ID.to.Bpifa2.tfs)
toptfs <- rbind((ID.to.Bpifa2.tfs %>% top_n(5,avg_logFC)), (ID.to.Bpifa2.tfs %>% top_n(-5,avg_logFC)))


#### Intercalated duct to/from SD 
ID.to.SD <- FindMarkers(object = epithelium, ident.1 = c("Adult_Intercalated duct"), ident.2 = "Adult_Striated duct", only.pos = F)
write.csv(ID.to.SD, file = "Transdiff potential - Intercalated duct to Striated duct.csv", row.names = T)
ID.to.SD.tfs <- ID.to.SD[rownames(ID.to.SD) %in% Transcription.factors$Gene, ]
write.csv(ID.to.SD.tfs, file = "Transdiff potential - Intercalated duct to Striated duct - Transcription Factors.csv", row.names = T)
ID.to.SD$gene <- rownames(ID.to.SD)
ID.to.SD.tfs$gene <- rownames(ID.to.SD.tfs)
toptfs <- rbind((ID.to.SD.tfs %>% top_n(5,avg_logFC)), (ID.to.SD.tfs %>% top_n(-5,avg_logFC)))

#### Intercalated duct to/from GCT 
ID.to.GCT <- FindMarkers(object = epithelium, ident.1 = c("Adult_Intercalated duct"), ident.2 = "Adult_GCT", only.pos = F)
write.csv(ID.to.GCT, file = "Transdiff potential - Intercalated duct to GCT.csv", row.names = T)
ID.to.GCT.tfs <- ID.to.GCT[rownames(ID.to.GCT) %in% Transcription.factors$Gene, ]
write.csv(ID.to.GCT.tfs, file = "Transdiff potential - Intercalated duct to GCT - Transcription Factors.csv", row.names = T)
ID.to.GCT$gene <- rownames(ID.to.GCT)
ID.to.GCT.tfs$gene <- rownames(ID.to.GCT.tfs)
toptfs <- rbind((ID.to.GCT.tfs %>% top_n(5,avg_logFC)), (ID.to.GCT.tfs %>% top_n(-5,avg_logFC)))

#### Intercalated duct to/from BD 
ID.to.BD <- FindMarkers(object = epithelium, ident.1 = c("Adult_Intercalated duct"), ident.2 = "Adult_Basal duct", only.pos = F)
write.csv(ID.to.BD, file = "Transdiff potential - Intercalated duct to BD.csv", row.names = T)
ID.to.BD.tfs <- ID.to.BD[rownames(ID.to.BD) %in% Transcription.factors$Gene, ]
write.csv(ID.to.BD.tfs, file = "Transdiff potential - Intercalated duct to BD - Transcription Factors.csv", row.names = T)
ID.to.BD$gene <- rownames(ID.to.BD)
ID.to.BD.tfs$gene <- rownames(ID.to.BD.tfs)
toptfs <- rbind((ID.to.BD.tfs %>% top_n(5,avg_logFC)), (ID.to.BD.tfs %>% top_n(-5,avg_logFC)))

#### Intercalated duct to/from Ascl3 
ID.to.Ascl3 <- FindMarkers(object = epithelium, ident.1 = c("Adult_Intercalated duct"), ident.2 = "Adult_Ascl3", only.pos = F)
write.csv(ID.to.Ascl3, file = "Transdiff potential - Intercalated duct to Ascl3.csv", row.names = T)
ID.to.Ascl3.tfs <- ID.to.Ascl3[rownames(ID.to.Ascl3) %in% Transcription.factors$Gene, ]
write.csv(ID.to.Ascl3.tfs, file = "Transdiff potential - Intercalated duct to Ascl3 - Transcription Factors.csv", row.names = T)
ID.to.Ascl3$gene <- rownames(ID.to.Ascl3)
ID.to.Ascl3.tfs$gene <- rownames(ID.to.Ascl3.tfs)
toptfs <- rbind((ID.to.Ascl3.tfs %>% top_n(5,avg_logFC)), (ID.to.Ascl3.tfs %>% top_n(-5,avg_logFC)))

######################################################################
######################################################################
################# CORRELATION ANALYSIS #########################
######################################################################
######################################################################

library(corrplot)
library(RColorBrewer)

matrix<-GetAssayData(acinar.epithelium)
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Bhlha15",])
correlations<-apply(matrix_mod,1,function(x){cor.test(gene,x)})
df <- data.frame(matrix(unlist(correlations), nrow=23737, byrow=T),stringsAsFactors=FALSE)
rownames(df) <- names(correlations)
colnames(df) <- names(correlations$Xkr4)
names(df)[3] <-"pvalue"
names(df)[4] <-"R"
df$pvalue <- as.numeric(df$pvalue)
df$R <- as.numeric(df$R)

write.csv(df, file = "Mist1-correlated genes.csv")
dftfs <- df[rownames(df) %in% Transcription.factors$Gene, ]

correlations2 <- as.data.frame(correlations)
correlations2$gene <- rownames(correlations2)
correlations3 <- correlations2[correlations2$gene %in% Transcription.factors$Gene, ]
correlations3 <- rbind(correlations3 %>% top_n(10,correlations), correlations3 %>% top_n(-10, correlations))
correlations3 <- correlations3[order(correlations3$correlations), ]

test <- as.matrix(GetAssayData(acinar.epithelium))
test <- subset.matrix(test, rownames(test) %in% Transcription.factors$Gene)
test <- as.matrix(t(test))
test.run <- cor(test)
test.run.subset <- subset.matrix(test.run, rownames(test.run) %in% correlations3$gene)
test.run.subset <- subset.matrix(test.run.subset, select = rownames(test.run) %in% correlations3$gene)

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
p.mat <- cor.mtest(test)
head(p.mat[, 1:5])
pmat.subset <- subset.matrix(p.mat, rownames(p.mat) %in% correlations3$gene)
pmat.subset <- subset.matrix(pmat.subset, select = rownames(p.mat) %in% correlations3$gene)


corrplot(test.run.subset, method="circle",tl.col = "black", order = "hclust",
         col = rev(brewer.pal(n=8, name="Spectral")),type = "upper",sig.level = 0.01, p.mat = pmat.subset, 
         insig = "blank") + theme(panel.background = element_blank())

corrplot(test.run.subset, method="color",bg = "transparent",diag = T,outline = "black",addgrid.col = "black",
         tl.col = "black", order = "FPC", 
         col = rev(brewer.pal(n=8, name="Spectral")),type = "upper",sig.level = 0.01, insig="pch", p.mat = pmat.subset) 

library(heatmaply)
heatmaply_cor(test.run.subset)

##### repeat for bpifa2 population
gene<-as.numeric(matrix_mod["Bpifa2",])
corr.bpifa2<-apply(matrix_mod,1,function(x){cor(gene,x)})
corr.bpifa2<- as.data.frame(corr.bpifa2)
corr.bpifa2$gene <- rownames(corr.bpifa2)
corr.bpifa2.tfs <- corr.bpifa2[corr.bpifa2$gene %in% Transcription.factors$Gene, ]
corr.bpifa2.tfs <- rbind(corr.bpifa2.tfs %>% top_n(10,corr.bpifa2), corr.bpifa2.tfs %>% top_n(-10, corr.bpifa2))

corr.bpifa2.tfs.matrix <- subset.matrix(test.run, rownames(test.run) %in% corr.bpifa2.tfs$gene)
corr.bpifa2.tfs.matrix <- subset.matrix(corr.bpifa2.tfs.matrix, select = rownames(test.run) %in% corr.bpifa2.tfs$gene)
pmat.bpifa.subset <- subset.matrix(p.mat, rownames(p.mat) %in% corr.bpifa2.tfs$gene)
pmat.bpifa.subset <- subset.matrix(pmat.bpifa.subset, select = rownames(p.mat) %in% corr.bpifa2.tfs$gene)

corrplot(corr.bpifa2.tfs.matrix, method="color",bg = "transparent",diag = T,outline = "black",addgrid.col = "black",
         tl.col = "black", order = "FPC", 
         col = brewer.pal(n=8, name="Spectral"),type = "upper",sig.level = 0.01, p.mat = pmat.subset, 
         insig = "pch") 


epithelium <- SetIdent(epithelium, value = "Epi.CellType")

plot1<-FeatureScatter(acinar.epithelium,feature1 = "Bhlha15",feature2 = "Aqp5", group.by = "Epi.CellType", pt.size = 1, cols = c9) +NoLegend()
plot1<-FeatureScatter(acinar.epithelium,feature1 = "Bhlha15",feature2 = "Ybx1", group.by = "Epi.CellType", pt.size = 1, cols = c9) +NoLegend()
plot1<-FeatureScatter(acinar.epithelium,feature1 = "Bhlha15",feature2 = "Ybx1", group.by = "Epi.CellType", pt.size = 1, cols = c9) +NoLegend()
plot1<-FeatureScatter(acinar.epithelium,feature1 = "Bhlha15",feature2 = "Ybx1", group.by = "Epi.CellType", pt.size = 1, cols = c9) +NoLegend()
plot1<-FeatureScatter(acinar.epithelium,feature1 = "Bhlha15",feature2 = "Sfrp1", group.by = "Epi.CellType", pt.size = 1, cols = c9) +NoLegend()

acinar.epithelium <- SetIdent(acinar.epithelium, value = "Epi.CellType.Stage")
DoHeatmap(acinar.epithelium, features = correlations3$gene, angle = 90,group.colors = c9) + scale_fill_gradientn(colors = c("navy", "darkseagreen", "yellow", "firebrick2", "firebrick4"))









