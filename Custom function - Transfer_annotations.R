
## create custom function to transfer annotations from one seurat object to another with same barcodes or cell identifiers
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
  
