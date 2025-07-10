library(data.table)
library(Matrix)

# Defined function
read_mat<-function(mat){
  data<-fread(mat,header=T)
  data$cellID<-paste(data$CellID)
  gene=unique(data$geneID)
  cell=unique(data$cellID)
  gene_idx=c(1:length(gene))#gene index   
  cell_idx=c(1:length(cell))#cell index
  names(gene_idx)=gene
  names(cell_idx)=cell
  mat=sparseMatrix(i=gene_idx[data$geneID],j=cell_idx[data$cellID],x=data$MIDCount)
  rownames(mat)=gene
  
  x=data$x[match(unique(data$CellID),data$CellID)]
  y=data$y[match(unique(data$CellID),data$CellID)]
  colnames(mat)=paste0(x,"_",y)
  return(mat)
}
load_obj <- function(mat){
  SeuObj<-CreateSeuratObject(counts = mat)
  SeuObj[["percent.mt"]] <- PercentageFeatureSet(SeuObj, pattern = "^mt-|^MT-|^Mt-")
  SeuObj@meta.data$coor_x=sub(rownames(SeuObj@meta.data),pattern = "_.*",replacement = "")
  SeuObj@meta.data$coor_y=sub(rownames(SeuObj@meta.data),pattern = ".*_",replacement = "")
  SeuObj@meta.data$coor_x=sub(SeuObj@meta.data$coor_x,pattern = "X",replacement = "")
  SeuObj@meta.data$coor_x=as.integer(SeuObj@meta.data$coor_x)
  SeuObj@meta.data$coor_y=as.integer(SeuObj@meta.data$coor_y)
  
  SeuObj
}

# Load raw data
raw_data =  'data/amy_spatial_raw/A02778D6.cellbin.gem'
mat <- read_mat(read_name)

# Create Seurat Object
seuobj <- load_obj(mat)

# Add spatial coords
loc = as.matrix(seuobj@meta.data[,c('coor_x','coor_y')])
colnames(loc) <- c('spatial_1','spatial_2')
loc = CreateDimReducObject(embeddings = loc)
seuobj@reductions$spatial <- loc

# Save data
saveRDS(seuobj,'data/amy_spatial/st_GW23_raw.rds')