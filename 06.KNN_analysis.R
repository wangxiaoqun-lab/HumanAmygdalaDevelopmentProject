library(Seurat)
library(dplyr)
library(htmlwidgets)
library(networkD3)

# read data
ge = readRDS('./DATA/ge_rename.rds')
InN = readRDS('./DATA/InN.rds')

# normalized data
temp.list <- list(ge<- ge,InN <- InN)
temp.list <- lapply(X = temp.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# get integrating features and anchors
features <- SelectIntegrationFeatures(object.list = temp.list)
anchors.InN <- FindIntegrationAnchors(object.list = temp.list, anchor.features = features,dims = 1:50)

InN.combined <- IntegrateData(anchorset = anchors.InN)
DefaultAssay(InN.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
InN.combined <- ScaleData(InN.combined, verbose = FALSE)
InN.combined <- RunPCA(InN.combined, npcs = 30, verbose = FALSE)

InN.combined <- RunPCA(InN.combined, npcs = 30, verbose = FALSE)
InN.combined <- RunUMAP(InN.combined, reduction = "pca", dims = 1:30,verbose = F,min.dist = 0.5)

InN.combined@meta.data$anno = paste0('InN ',sprintf("%02s", InN.combined@meta.data$leiden))
InN.combined@meta.data$anno[InN.combined@meta.data$anno == 'InN NA'] = InN.combined@meta.data$celltype_new3[InN.combined@meta.data$anno == 'InN NA']
Idents(InN.combined) <- 'anno'

options(repr.plot.width=12,repr.plot.height=10)
DimPlot(InN.combined,label = T,pt.size = 2)

# KNN prediction
InN.combined   <-   FindNeighbors(InN.combined,k.param = 15, reduction = "pca", dims = 1:30)
snn <- InN.combined@graphs$integrated_snn
labels <-InN.combined$anno
n <- nrow(snn)
pred_label <- rep(NA, n)
pred_conf  <- rep(NA, n)
for (i in 1:n) {
  neigh <- which(snn[i, ] > 0)
  neigh <- setdiff(neigh, i)
  if (length(neigh) == 0) next
  weights <- snn[i, neigh]
  neigh_labels <- labels[neigh]
  score <- tapply(weights, neigh_labels, sum)
  score_sorted <- sort(score, decreasing = TRUE)
  prob <- score / sum(score)
  if (length(score_sorted) > 1 && score_sorted[1] == score_sorted[2]) {
    pred_label[i] <- NA
    pred_conf[i] <- max(prob)
    next
  }
  pred_label[i] <- names(score_sorted)[1]
  pred_conf[i]  <- max(prob)
}

# save result
InN.combined$predicted_celltype <- pred_label

# prepare sankey plot
plot_df = InN.combined@meta.data[InN.combined@meta.data$orig.ident=='AMY',]
plot_df = plot_df[! is.na(plot_df$predicted_celltype),]
GE_label = c('Progenitor','LGE','CGE','MGE')
plot_df_sub1 = plot_df[plot_df$predicted_celltype2 %in% GE_label,]
flow_df1 <- plot_df_sub1 %>%
  count(Branch, predicted_celltype) %>%
  group_by(Branch) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()
plot_df_sub2 = plot_df[! plot_df$predicted_celltype2 %in% GE_label,]
flow_df2 <- plot_df_sub2 %>%
  count(Branch, predicted_celltype) %>%
  group_by(Branch) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

# merge df
flow_df = rbind(flow_df1,flow_df2)
nodes <- data.frame(
  name = unique(c(flow_df$predicted_celltype, flow_df$Branch))
)

# get links
links <- flow_df %>%
  mutate(
    source = match(predicted_celltype, nodes$name) - 1,
    target = match(Branch, nodes$name) - 1,
    value = percent   
  )

# plot
fig = sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  fontSize = 12,
  nodeWidth = 30,
 width = 600,   
  height = 600
)