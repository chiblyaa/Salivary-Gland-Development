
library(Seurat)
## This script was only used to transfer final annotations to the integrated files

#load integrated files
smg.p.integrated <- readRDS("~/10X Manuscript Projects/10X revisions/Postnatal SMG Integrated (P1-Adult).rds")
smg.e.integrated <- readRDS("~/10X Manuscript Projects/10X revisions/Embryonic SMG Integrated (E12-E16).rds")

##load files with fixed annotations and updated metadata
e12epi <- readRDS("~/10X Manuscript Projects/10X revisions/E12-E16 epithelium/E12 epithelium qcd_annotated.rds")
e14epi <- readRDS("~/10X Manuscript Projects/10X revisions/E12-E16 epithelium/E14 epithelium qcd_annotated.rds")
e16epi <- readRDS("~/10X Manuscript Projects/10X revisions/E12-E16 epithelium/E16 epithelium qcd_annotated.rds")
p1epi <- readRDS("~/10X Manuscript Projects/10X revisions/P1-Adult/p1 epithelium qcd_annotated.rds")
p30epi <- readRDS("~/10X Manuscript Projects/10X revisions/P1-Adult/p30 epithelium qcd_annotated.rds")
adultepi <- readRDS("~/10X Manuscript Projects/10X revisions/P1-Adult/adult epithelium qcd_annotated.rds")

## create custom function to transfer annotations
transfer.annotation <- function(to, from, metadata){
  require(Seurat)
  for (j in 1:length(metadata)) {
    for (i in 1:length(from)) {
      from[[i]] <- SetIdent(from[[i]], value = metadata[j])
      from.cells <- WhichCells(from[[i]])
      Idents(object = to, cells=from.cells) <- Idents(from[[i]])
    }
    to[[metadata[j]]] <- Idents(to)
  }
  return(to)
}
  
smg.e.integrated <- transfer.annotation(to = smg.e.integrated, from = list(e12epi, e14epi, e16epi), metadata = c("celltype.fixed", "subpopulations"))
smg.p.integrated <- transfer.annotation(to = smg.p.integrated, from = list(p1epi, p30epi, adultepi), metadata = c("celltype.fixed", "subpopulations"))

cell.type.levels <- c("End bud", "Krt19+ duct", "Basal duct", "Myoepithelial", "Bpifa2+ Proacinar", "Smgc+ Proacinar", "Mitotic cells",
                      "Serous acinar", "Seromucous acinar", "Gstt1+ Female", "Gstt1+ Male", "Gstt1+", "Intercalated duct", "Ascl3+ duct", "GCT", "Striated duct",    
                      "Endothelial", "Smooth muscle", "Nerves",  "Glial cells",   "Macrophages", "Mast cells", "NK cells","Mesenchyme", "Stromal", "Erythroid", "Undefined") # use this for consistent sorting of cell types in plots (optional)

# Rename identities based on our findings
smg.p.integrated <- SetIdent(smg.p.integrated, value = "celltype.fixed")
smg.p.integrated <- RenameIdents(smg.p.integrated, "Smgc+ Female" = "Gstt1+ Female", "Smgc+" = "Gstt1+ Female", 
                                 "Smgc+ Male" = "Gstt1+ Male", "Acinar" = "Seromucous acinar",
                                 "Bpifa2+" = "Serous acinar")
DimPlot(smg.p.integrated)  
Idents(smg.p.integrated) <- factor(Idents(smg.p.integrated), levels = cell.type.levels)
smg.p.integrated$celltype.fixed <- Idents(smg.p.integrated)
DimPlot(smg.p.integrated, group.by = "celltype.fixed")
saveRDS(smg.p.integrated, file = "Postnatal SMG Integrated - final annotations.rds")

smg.e.integrated <- SetIdent(smg.e.integrated, value = "celltype.fixed")
DimPlot(smg.e.integrated)  
Idents(smg.e.integrated) <- factor(Idents(smg.e.integrated), levels = cell.type.levels)
smg.e.integrated$celltype.fixed <- Idents(smg.e.integrated)
DimPlot(smg.e.integrated, group.by = "celltype.fixed")  
saveRDS(smg.e.integrated, file = "Embryonic SMG Integrated - final annotations.rds")

#### Split individual objects with final annotations for SGMap
embryonic.list <- SplitObject(smg.e.integrated, split.by = "stage")
postnatal.list <- SplitObject(smg.p.integrated, split.by = "stage")

### Generate list of cell-type genes with final annotations for each stage

cell.DEGs <- list()
for (i in 1:length(embryonic.list)) {
  DefaultAssay(embryonic.list[[i]]) <- "RNA"
  embryonic.list[[i]] <- SetIdent(embryonic.list[[i]], value= "celltype.fixed")
  cell.DEGs[[i]] <- FindAllMarkers(embryonic.list[[i]], only.pos = T, min.pct = 0.2)  
}

cell.DEGs.postnatal <- list()
for (i in 1:length(postnatal.list)) {
  DefaultAssay(postnatal.list[[i]]) <- "RNA"
  postnatal.list[[i]] <- SetIdent(postnatal.list[[i]], value= "celltype.fixed")
  cell.DEGs.postnatal[[i]] <- FindAllMarkers(postnatal.list[[i]], only.pos = T, min.pct = 0.2)  
}

write.csv(x = cell.DEGs[[1]], file = "E12 cell type markers (final annotations).csv")
write.csv(x = cell.DEGs[[2]], file = "E14 cell type markers (final annotations).csv")
write.csv(x = cell.DEGs[[3]], file = "E16 cell type markers (final annotations).csv")
write.csv(x = cell.DEGs.postnatal[[1]], file = "P1 cell type markers (final annotations).csv")
write.csv(x = cell.DEGs.postnatal[[2]], file = "P30 cell type markers (final annotations).csv")
write.csv(x = cell.DEGs.postnatal[[3]], file = "Adult cell type markers (final annotations).csv")

### generate figures
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
            "hotpink",
            "#0074c9", 
            "#7f4a21", 
            "#c3c1ff",
            "#fbe54b")

cell.types <- c("Seromucous acinar", "Ascl3+ duct", "Basal duct", "Serous acinar", "Bpifa2+ Proacinar", "End bud", "Endothelial",
                "Erythroid", "GCT", "Glial cells", "Intercalated duct", "Krt19+ duct", "Macrophages", "Mast cells",
                "Mesenchyme", "Mitotic cells", "Myoepithelial", "Nerves", "NK cells", "Gstt1+ Female", "Gstt1+ Male",
                "Smgc+ Proacinar", "Smooth muscle", "Striated duct", "Stromal")

names(colors) <- cell.types

pdf("UMAPs with final annotations.pdf", useDingbats = F, width = 6.5, height = 5)
DimPlot(smg.p.integrated, group.by = "celltype.fixed", cols = colors)
DimPlot(smg.e.integrated, group.by = "celltype.fixed", cols = colors)
dev.off()
