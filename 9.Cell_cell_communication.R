library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# Load data
all <- readRDS('data/AMY/result.rds')

# Prepare data for CellChat
temp = subset(all,cellchat_sub == TRUE)
data.input <- temp@assays$RNA@data
meta = temp@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))

# Compute cell-cell communication net
cellchat@DB <- CellChatDB.human
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat,raw.use=FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Plot
print(netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions"))
print(netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength"))
netVisual_chord_gene(cellchat, sources.use = c(1,3,4), targets.use = c(2), lab.cex = 0.5,legend.pos.y = 30)