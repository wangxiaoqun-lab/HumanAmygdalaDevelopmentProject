library(Seurat)
library(Matrix)
library(FNN)
library(reshape2)
library(scales)
library(ggplot2)
library(dplyr)
library(gplots)
source("./denpendency/help_code.R")

# Load data
InN <- readRDS('data/AMY/result.rds')

# GE to InN
time_point = c ('GE','GW8-10')
time_i = time_point[1]
time_j = time_point[2]
print(time_i)
print(time_j)
emb = readRDS(paste0("data/AMY/integrated_", time_i, "_", time_j, "_umap3.rds"))
emb = data.frame(emb)
anno1 <- readRDS('data/GE.rds')
anno1$anno <- paste0('GE:',ge_data$subtype_new)
anno1$Anno = as.vector(anno1$anno)
anno1$day = "pre"
anno1 = anno1[[]][,c("day", "Anno")]
rownames(anno1) <- paste0(rownames(anno1),'_1')
anno2 = subset(InN,stage==time_j)
anno2$anno5 <- paste0(anno2$stage,':',anno2$anno4)
anno2$Anno = as.vector(anno2$anno5)
anno2$day = "nex"
anno2 = anno2[[]][,c("day", "Anno")]
rownames(anno2) <- paste0(rownames(anno2),'_2')
anno = rbind(anno1, anno2)

if(nrow(emb) != nrow(anno)){
        print("Error!")
    }
anno = anno[rownames(emb),]
res = createLineage_Knn(emb, anno,  k_neigh = 5)
saveRDS(res, paste0("data/AMY/integrated_", time_i, "_", time_j, "_Knn_umap.rds"))
for(k_i in c(10, 20, 30)){
    res = createLineage_Knn(emb, anno,  k_neigh = k_i) 
    saveRDS(res, paste0("data/AMY/integrated_", time_i, "_", time_j, "_Knn_umap_k_", k_i, ".rds"))
}
replication_times=500
dat = res
state_1 = row.names(dat[[1]])
state_2 = names(dat[[1]])
tmp_1 = matrix(NA,nrow(dat[[1]]),ncol(dat[[1]]))
for(i in 1:nrow(dat[[1]])){
    for(j in 1:ncol(dat[[1]])){
        xx = NULL
        for(k in 1:replication_times){
            xx = c(xx, dat[[k]][i,j])
        }
    tmp_1[i,j] = median(xx[!is.na(xx)])
    }
}
tmp_1 = data.frame(tmp_1)
row.names(tmp_1) = state_1
names(tmp_1) = state_2
write.csv(tmp_1, paste0("./integrated_", time_i, "_", time_j, "_Knn_umap.csv"))

# Fetal InN
time_point=c('GW8-10','GW11-14','GW17-21','GW22-25')
for (kk in c(1,2,3)){
    time_i = time_point[kk]
    time_j = time_point[kk+1]
    print(time_i)
    print(time_j)
    emb = readRDS(paste0("data/AMY/integrated_", time_i, "_", time_j, "_umap3.rds"))
    emb = data.frame(emb)
    anno1 = subset(InN,stage==time_i)
    anno1$anno5 <- paste0(anno1$stage,':',anno1$anno4)
    anno1$Anno = as.vector(anno1$anno5)
    anno1$day = "pre"
    anno1 = anno1[[]][,c("day", "Anno")]
    anno2 = subset(InN,stage==time_j)
    anno2$anno5 <- paste0(anno2$stage,':',anno2$anno4)
    anno2$Anno = as.vector(anno2$anno5)
    anno2$day = "nex"
    anno2 = anno2[[]][,c("day", "Anno")]
    anno = rbind(anno1, anno2)
    if(nrow(emb) != nrow(anno)){
        print("Error!")
    }
    anno = anno[rownames(emb),]
    res = createLineage_Knn(emb, anno,  k_neigh = 5)
    saveRDS(res, paste0("data/AMY/integrated_", time_i, "_", time_j, "_Knn_umap.rds"))
    for(k_i in c(10, 20, 30)){
        res = createLineage_Knn(emb, anno,  k_neigh = k_i) #### createLineage_Knn function was in help_code.R
        saveRDS(res, paste0("data/AMY/integrated_", time_i, "_", time_j, "_Knn_umap_k_", k_i, ".rds"))
    }
    replication_times=500
    dat = res
    state_1 = row.names(dat[[1]])
    state_2 = names(dat[[1]])
    
    tmp_1 = matrix(NA,nrow(dat[[1]]),ncol(dat[[1]]))
    for(i in 1:nrow(dat[[1]])){
        for(j in 1:ncol(dat[[1]])){
            xx = NULL
            for(k in 1:replication_times){
                xx = c(xx, dat[[k]][i,j])
            }
        tmp_1[i,j] = median(xx[!is.na(xx)])
        }
    }
    tmp_1 = data.frame(tmp_1)
    row.names(tmp_1) = state_1
    names(tmp_1) = state_2

    write.csv(tmp_1, paste0("./integrated_", time_i, "_", time_j, "_Knn_umap.csv"))
}

# Fetal to adult
time_point = c ('GW22-25','Adult')
time_i = time_point[1]
time_j = time_point[2]
print(time_i)
print(time_j)
emb = readRDS(paste0("data/AMY/integrated_", time_i, "_", time_j, "_umap3.rds"))
emb = data.frame(emb)
anno1 = subset(InN,stage==time_i)
anno1$anno5 <- paste0(anno1$stage,':',anno1$anno4)
anno1$Anno = as.vector(anno1$anno5)
anno1$day = "pre"
anno1 = anno1[[]][,c("day", "Anno")]
anno2 <- subset(Hu_orig_onth,celltype=='InN')
anno2$anno <- gsub(pattern = 'Human_',replacement = '',x = anno2$orig_anno)
anno2$anno2 <- paste0('Adult:',anno2$anno)
anno2$Anno = as.vector(anno2$anno2)
anno2$day = "nex"
anno2 = anno2[[]][,c("day", "Anno")]
anno = rbind(anno1, anno2)
if(nrow(emb) != nrow(anno)){
        print("Error!")
    }
anno = anno[rownames(emb),]
res = createLineage_Knn(emb, anno,  k_neigh = 5)
saveRDS(res, paste0("data/AMY/integrated_", time_i, "_", time_j, "_Knn_umap.rds"))
for(k_i in c(10, 20, 30)){
    res = createLineage_Knn(emb, anno,  k_neigh = k_i) #### createLineage_Knn function was in help_code.R
    saveRDS(res, paste0("data/AMY/integrated_", time_i, "_", time_j, "_Knn_umap_k_", k_i, ".rds"))
}
replication_times=500
dat = res
state_1 = row.names(dat[[1]])
state_2 = names(dat[[1]])
tmp_1 = matrix(NA,nrow(dat[[1]]),ncol(dat[[1]]))
for(i in 1:nrow(dat[[1]])){
    for(j in 1:ncol(dat[[1]])){
        xx = NULL
        for(k in 1:replication_times){
            xx = c(xx, dat[[k]][i,j])
        }
    tmp_1[i,j] = median(xx[!is.na(xx)])
    }
}
tmp_1 = data.frame(tmp_1)
row.names(tmp_1) = state_1
names(tmp_1) = state_2
write.csv(tmp_1, paste0("./integrated_", time_i, "_", time_j, "_Knn_umap.csv"))

# Summary
replication_times=500
res_median_umap = list()
time_point = c('GE','GW8-10','GW11-14','GW17-21','GW22-25','Adult')
for(time_i in 1:(length(time_point)-1)){
  print(paste0(time_point[time_i], ":", time_point[time_i+1]))
  dat = readRDS(paste0("data/AMY/integrated_",time_point[time_i],"_",time_point[time_i+1],"_Knn_umap.rds"))
    
  state_1 = row.names(dat[[1]])
  state_2 = names(dat[[1]])
  tmp_1 = matrix(NA,nrow(dat[[1]]),ncol(dat[[1]]))
  for(i in 1:nrow(dat[[1]])){
    for(j in 1:ncol(dat[[1]])){
      xx = NULL
      for(k in 1:replication_times){
        xx = c(xx, dat[[k]][i,j])
      }
      tmp_1[i,j] = median(xx[!is.na(xx)])
    }
  }
  tmp_1 = data.frame(tmp_1)
  row.names(tmp_1) = state_1
  names(tmp_1) = state_2
  res_median_umap[[time_i]] = tmp_1
}

dat = NULL
for(i in 1:length(res_median_umap)){
  print(time_point[i])
  dat = rbind(dat, melt(as.matrix(res_median_umap[[i]])))
}

dat = data.frame(dat)
names(dat) = c("nex", "pre", "prob")
dat$pre_time = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][1]))
dat$pre_cell = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][2]))
dat$nex_time = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][1]))
dat$nex_cell = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][2]))

# Filter data
InN$group <- paste0(InN$stage,':',InN$anno4)
remove_gourp <- names(table(InN$group))[table(InN$group) < 100]
dat <- dat[(!dat$nex %in% remove_gourp) & (!dat$pre %in% remove_gourp),]
dat$nex <- as.character(dat$nex)
dat$pre <- as.character(dat$pre)
dat$pre <- plyr::mapvalues(x = dat$pre,from = 'GE:Progenitor',to = 'GE_pre:Progenitor')

# Add root
dat$pre_time = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][1]))
dat$pre_cell = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][2]))
dat$nex_time = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][1]))
dat$nex_cell = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][2]))
saveRDS(dat, paste0("data/AMY/integrated_edge_all.rds"))
print(paste0("how many edges: ", nrow(dat)))
print(paste0("how many edges (> 0): ", nrow(dat[dat$prob>0,])))
print(paste0("how many edges (> 0.1): ", nrow(dat[dat$prob>=0.1,])))
print(paste0("how many edges (> 0.2): ", nrow(dat[dat$prob>=0.2,])))
print(paste0("how many edges (> 0.7): ", nrow(dat[dat$prob>=0.7,])))
print(paste0("how many edges (> 0.8): ", nrow(dat[dat$prob>=0.8,])))
print(paste0("how many nodes: ", length(unique(c(as.vector(dat$pre), as.vector(dat$nex))))))
print(paste0("how many cell types: ", length(unique(c(as.vector(dat$pre_cell), as.vector(dat$nex_cell))))))
x = dat[dat$prob>=0.1,]
x = x[,c("pre","nex","prob")]
print(paste0("how many nodes now: ", length(unique(c(as.vector(x$pre), as.vector(x$nex))))))
dummy = NULL
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:C1',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:C2',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:C3',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:L1',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:L2',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:L3',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:L4',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:L5',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:L6',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:L7',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:M1',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:M2',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:M3',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:M4',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:M5',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:M6',1))
dummy = rbind(dummy, c("GE_pre:Progenitor",'GE:M7',1))

dummy = data.frame(dummy)
names(dummy) = c("pre","nex","prob")

res = rbind(x, dummy)

res$pre <- as.character(res$pre)
res$nex <- as.character(res$nex)
res$prob <- as.double(res$prob)

dat_sub = res
dat_sub$pre_cell = unlist(lapply(as.vector(dat_sub$pre), function(x) strsplit(x,"[:]")[[1]][2]))
dat_sub$nex_cell = unlist(lapply(as.vector(dat_sub$nex), function(x) strsplit(x,"[:]")[[1]][2]))
sum(dat_sub$pre_cell == dat_sub$nex_cell)
sum(dat_sub$pre_cell == dat_sub$nex_cell)/nrow(dat_sub)
sum(dat_sub$pre_cell != dat_sub$nex_cell)
sum(dat_sub$pre_cell != dat_sub$nex_cell)/nrow(dat_sub)

# Statistic
print(paste0("how many edges: ", nrow(res)))
print(paste0("how many nodes: ", length(unique(c(as.vector(res$pre), as.vector(res$nex))))))

nex_list = as.vector(unique(res$nex))
tree = NULL
for(i in 1:length(nex_list)){
  res_sub = res[res$nex==nex_list[i],]
  if(nrow(res_sub)==1){
    tree = rbind(tree, res_sub)
  } else {
    res_sub = res_sub[order(res_sub$prob, decreasing = TRUE),]
    tree = rbind(tree, res_sub[1,])
  }
}
tree = data.frame(tree)
tree = tree[,c("pre","nex")]

write.table(res, paste0("./integrated_edge_prob.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(tree, paste0("./integrated_edge.txt"), row.names = F, col.names = F, quote = F, sep = "\t")