library(Seurat)

# read data
sc_data = readRDS('./DATA/sc_stage1.rds')
st_gw10 = readRDS('./DATA/st_gw10.rds')

# transfer label function
transfer_data <- function(data1,data2){

    # normalize data
    data1 <- FindVariableFeatures(data1)
    data1 <- ScaleData(data1)
    data2 <- NormalizeData(data2) 
    data2 <- FindVariableFeatures(data2)
    data2 <- ScaleData(data2)
	
	# get anchors
    transfer.anchors <- FindTransferAnchors(reference = data2, query = data1,dims = 1:30,
    reference.assay = "RNA", query.assay = "RNA", reduction = "cca")
    
	# get celltype predictions
    celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = data2$anno, weight.reduction = 'cca', dims = 1:30)
    data1 <- AddMetaData(data1, metadata = celltype.predictions['predicted.id'],col.name = 'anno')
    data1 <- AddMetaData(data1, metadata = celltype.predictions['prediction.score.max'],col.name = 'transfer_score')
    return (data1)
}

# save predicted annotation
st_gw10 = transfer_data(data1 = st_gw10,data2 = sc_data)
write.csv(st_gw10@meta.data,'./meta_st_10_slice1.csv')