library(Seurat)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(VGAM)
library(parallel)
library(circlize)
library(ComplexHeatmap)
library(dendsort)
library(ggplotify)
library(viridis)
library(RColorBrewer)

# read progenitor data
Pro <- readRDS('./DATA/Pro.rds')
Pro_sub <- subset(Pro,anno %in% c('RG1','RG2','RG3','RG4','RG5','ExN IP1','ExN IP2','InN IP1'))

# read whole brain data
adata <- readRDS('./DATA/human_dev_sub_tele_dien.rds')
adata <- NormalizeData(adata)
adata_sub  = subset(adata ,downsample=10000)

# get anchors
Pro_sub <- FindVariableFeatures(Pro_sub,nfeatures = 2000)
adata_sub <- FindVariableFeatures(adata_sub ,nfeatures = 2000)
anchors <- FindTransferAnchors(
  reference = adata_sub,
  query = Pro_sub,
  reduction = "cca",
  dims = 1:10
)

# transfer label
predictions <- TransferData(anchorset = anchors,weight.reduction = 'cca', refdata = adata_sub$Region,dims = 1:10)
Pro_sub <- AddMetaData(Pro_sub, metadata = predictions['predicted.id'],col.name = 'Region')

# visualization
mat = as.matrix(table(Pro_sub$anno,Pro_sub$Region))
mat =  mat / rowSums(mat)
options(repr.plot.width=3,repr.plot.height=6)
fig = pheatmap::pheatmap(mat[c('RG1','RG2','RG3','RG4','RG5','ExN IP1','ExN IP2','InN IP1'),c('Telencephalon','Diencephalon')]
                         ,cluster_cols = F,cluster_rows = F)
ggsave(filename = './FIGS/AMY_Pro_Cortex_heatmap.pdf',fig,width = 4,height = 6)

# RG4 lineage
RG4_sub = subset(Pro,anno %in% c('RG4','ExN IP2','Glia lineage2'))
Idents(RG4_sub) <- 'anno'
deg <- FindAllMarkers(RG4_sub)

# plot DEG
fig <- markerVocalno(markers = deg,topn = 20,log2FC = 0.5,labelCol = c('#ff7f0e', '#17becf', '#d62728'))+NoLegend()
ggsave('./FIGS/markerVolcano_RG4.pdf',fig,width = 10,height = 10)

# get monocle3 pseudotime
RG4_sub.monocle <- as.cell_data_set(RG4_sub)
RG4_sub.monocle <- cluster_cells(RG4_sub.monocle,resolution = 0.005)
plot_cells(RG4_sub.monocle,show_trajectory_graph = FALSE,color_cells_by = 'cluster',cell_size = 2)
RG4_sub.monocle <- learn_graph(RG4_sub.monocle,close_loop = FALSE,
  learn_graph_control = list(
    prune_graph = TRUE,
    rann.k = 30
  ))
RG4_sub.monocle <- order_cells(RG4_sub.monocle, root_pr_nodes=c('Y_76'))
RG4_sub$pseudotime <- pseudotime(RG4_sub.monocle)
plot_cells(RG4_sub.monocle, show_trajectory_graph = TRUE,trajectory_graph_segment_size = 2,color_cells_by = 'anno',,cell_size = 2)+
scale_color_manual(values = c(
  "ExN IP2" = '#ff7f0e', 
  "Glia lineage2" = '#d62728', 
  "RG4" = '#17becf'
))
ggsave('./FIGS/RG4_monocle.pdf',width=6,height = 6)

FeaturePlot(RG4_sub,features = 'pseudotime',pt.size=1)+scale_color_viridis_c(option = "plasma")
ggsave('./FIGS/RG4_pseudotime.pdf',width=6,height = 6)

# VGAM fittness
compareModels <- function(full_models, reduced_models){
  stopifnot(length(full_models) == length(reduced_models))
  test_res <- mapply(function(x,y) { 
    if (is.null(x) == FALSE && is.null(y) == FALSE) {
      lrt <- VGAM::lrtest(x,y) 
      pval=lrt@Body["Pr(>Chisq)"][2,]
      family = x@family@vfamily
      if (length(family) > 1)
        family = family[1]
      data.frame(status = "OK", family=family, pval=pval)
    } else { data.frame(status = "FAIL", family=NA, pval=1.0) } 
  } , full_models, reduced_models, SIMPLIFY=FALSE, USE.NAMES=TRUE)
  
  test_res <- do.call(rbind.data.frame, test_res)
  test_res$qval <- p.adjust(test_res$pval, method="BH")
  test_res
}
esponseMatrix <- function(models, newdata = NULL, response_type="response", cores = 1) {
  res_list <- mclapply(models, function(x) {
    if (is.null(x)) { NA } else {
      if (x@family@vfamily %in% c("negbinomial", "negbinomial.size")) {
        predict(x, newdata = newdata, type = response_type)
      } else if (x@family@vfamily %in% c("uninormal")) {
        predict(x, newdata = newdata, type = response_type)
      }
      else {
        10^predict(x, newdata = newdata, type = response_type)
      }
    }
  }, mc.cores = cores)
  
  res_list_lengths <- lapply(res_list[is.na(res_list) == FALSE],
                             length)
  stopifnot(length(unique(res_list_lengths)) == 1)
  num_na_fits <- length(res_list[is.na(res_list)])
  if (num_na_fits > 0) {
    na_matrix <- matrix(rep(rep(NA, res_list_lengths[[1]]),
                            num_na_fits), nrow = num_na_fits)
    row.names(na_matrix) <- names(res_list[is.na(res_list)])
    non_na_matrix <- Matrix::t(do.call(cbind, lapply(res_list[is.na(res_list) ==
                                                                FALSE], unlist)))
    row.names(non_na_matrix) <- names(res_list[is.na(res_list) ==
                                                 FALSE])
    res_matrix <- rbind(non_na_matrix, na_matrix)
    res_matrix <- res_matrix[names(res_list), ]
  }
  else {
    res_matrix <- Matrix::t(do.call(cbind, lapply(res_list, unlist)))
    row.names(res_matrix) <- names(res_list[is.na(res_list) ==
                                              FALSE])
  }
  res_matrix
}

# get top50 DEG
deg_top50 <- deg %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 100)
glia_genes = deg_top50[deg_top50$cluster == 'Glia lineage2',]$gene
ExN_genes = deg_top50[deg_top50$cluster == 'ExN IP2',]$gene
gene = c(glia_genes,ExN_genes)

# fittness APC part
temp <- subset(RG4_sub,anno=='Glia lineage2')
newdata1 <- data.frame(pseudotime = seq(min(temp$pseudotime),max(temp$pseudotime),length.out=100))
exp_mt <- as.matrix(temp@assays$RNA@scale.data)
cds_data <- data.frame(pseudotime=as.numeric(temp$pseudotime),as.data.frame(t(exp_mt[gene,])))
gene <- gsub(pattern = '-','_',gene)
colnames(cds_data) <- gsub(pattern = '\\.',replacement = '_',x = colnames(cds_data))

forma_list=list()
for (i in gene){
    forma <- paste0(i,'~sm.ns(pseudotime, df=3)')
    forma_list[[i]] <- vglm(formula =forma ,family = uninormal,data = cds_data, trace = TRUE)
}

smooth_mat_APC = esponseMatrix(models = forma_list,newdata = newdata1)
rownames(smooth_mat_APC) <- gsub('_','-',rownames(smooth_mat_APC))
col_anno <- columnAnnotation(Pseudotime=as.numeric(colnames(smooth_mat_APC)),col=list(Pseudotime=colorRamp2(c(0,100), c("white", "green"))))
row_anno <- rowAnnotation(DEG = c(rep('Glia',times = 50),rep('ExN',times = 50)))
temp1=dendsort(hclust(dist(smooth_mat_APC[1:50,])))

# fittness ExN part
temp <- subset(RG4_sub,anno=='ExN IP2')
gene = c(glia_genes,ExN_genes)
exp_mt <- as.matrix(temp@assays$RNA@scale.data)
cds_data <- data.frame(pseudotime=as.numeric(temp$pseudotime),as.data.frame(t(exp_mt[gene,])))
gene <- gsub(pattern = '-','_',gene)
colnames(cds_data) <- gsub(pattern = '\\.',replacement = '_',x = colnames(cds_data))

forma_list=list()
for (i in gene){
    forma <- paste0(i,'~sm.ns(pseudotime, df=3)')
    forma_list[[i]] <- vglm(formula =forma ,family = uninormal,data = cds_data, trace = TRUE)
}
smooth_mat_ExN = esponseMatrix(models = forma_list,newdata = newdata1)
rownames(smooth_mat_ExN) <- gsub('_','-',rownames(smooth_mat_ExN))
col_anno <- columnAnnotation(Pseudotime=as.numeric(colnames(smooth_mat_ExN)),col=list(Pseudotime=colorRamp2(c(0,100), c("white", "green"))))
row_anno <- rowAnnotation(DEG = c(rep('Glia',times = 50),rep('ExN',times = 50)))
temp2=dendsort(hclust(dist(smooth_mat_ExN[51:100,])))
fig <- Heatmap(smooth_mat_ExN,row_order=c(temp1$labels[temp1$order],temp2$labels[temp2$order]),cluster_rows = F,cluster_columns = F,show_column_names = F,left_annotation = row_anno,top_annotation = col_anno,use_raster = T)
fig <- as.ggplot(fig)

# merge matrix and plot
smooth_mat <- cbind(smooth_mat_ExN[,ncol(smooth_mat_ExN):1],smooth_mat_APC)
col2 <- circlize::colorRamp2(c(-1,-0.5,0,0.5,1),colors = brewer.pal(5, 'Spectral')[length(brewer.pal(5, 'Spectral')):1])
row_anno <- rowAnnotation('DEG' = factor(c(rep('APC',times = 50),rep('ExN',times = 50)),levels = c('ExN','APC')),col=list('DEG'=c('ExN'='#FF9858','APC'='#00CCEE')))
col_anno <- columnAnnotation(Pseudotime=as.numeric(colnames(smooth_mat)),col=list(Pseudotime=colorRamp2(seq(from = 0,to = 100,length.out = 20), plasma(20))))
fig <- Heatmap(smooth_mat,col = col2,cluster_rows = F,show_row_names = F,row_order=rev(c(temp1$labels[temp1$order],temp2$labels[temp2$order])),cluster_columns = F,show_column_names = F,left_annotation = row_anno,top_annotation = col_anno,use_raster = T,name = 'Gene expression')
fig <- as.ggplot(fig)
ggsave(filename = './FIGS/HEATMAP_RG4_lineage.pdf',plot = fig,width = 8,height = 8)