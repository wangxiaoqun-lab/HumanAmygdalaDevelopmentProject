library(SCENIC)
library(Seurat)
library(AUCell)
library(SCopeLoomR)

# Defined function
binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}

# Load data
loom <- open_loom('SCENIC_output/pyscenic_output.loom')

# Get regulon details
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name = 'RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
temp <- names(regulonAucThresholds)
names(temp) <- as.vector(regulonAucThresholds)
regulonAucThresholds <- temp
binaryRegulonActivity <- binarizeAUC(regulonAUC, regulonAucThresholds)

# Merge Seurat object
data[['regulon']] <- CreateAssayObject(counts = getAUC(regulonAUC))
data[['regulon_binary']] <- CreateAssayObject(counts = binaryRegulonActivity)
saveRDS(object = data,'data/AMY/SCENIC_result.rds')