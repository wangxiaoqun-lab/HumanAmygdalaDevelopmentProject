library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)

# Transfer data format
Convert("data/AMY/result.h5ad", dest = "h5seurat", overwrite = TRUE)
all <- LoadH5Seurat("data/AMY/result.h5seurat",assays = 'RNA')
saveRDS(all,'data/AMY/result.rds')

# Set category and Find DEGs
Idents(all) <- 'stage'
deg_stage <- FindAllMarkers(all,only.pos = T)

# Get gene id
shared_deg <- deg_stage[deg_stage$gene %in% VariableFeatures(all),]
temp <- bitr(shared_deg$gene,fromType = 'SYMBOL',toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")

# Gene enrichment
go_result <- compareCluster(ENTREZID~cluster,data=shared_deg,fun='enrichGO',OrgDb='org.Hs.eg.db',ont='BP')

# Set category
go_result@compareClusterResult$Cluster <- factor(go_result@compareClusterResult$Cluster,levels = c('GW8-10','GW11-14','GW17-21','GW22-25'))
go_result@compareClusterResult <- go_result@compareClusterResult[order(go_result@compareClusterResult$Cluster),]

# Choose GO terms
cate <- c('chromosome segregation','regulation of mitotic cell cycle','neural precursor cell proliferation',
         'GABAergic neuron differentiation','gliogenesis','glial cell development','oligodendrocyte development',
         'regulation of trans-synaptic signaling','axon extension','axon extension involved in axon guidance','learning or memory',
          'synapse organization','regulation of membrane potential','regulation of synaptic transmission, glutamatergic','excitatory postsynaptic potential','presynapse assembly','dendrite development'
         )

# Plot
fig = dotplot(object = go_result,showCategory=cate)
ggsave(filename = 'dotplot_GO_stage.pdf',width = 7,height = 8,dpi=400)
