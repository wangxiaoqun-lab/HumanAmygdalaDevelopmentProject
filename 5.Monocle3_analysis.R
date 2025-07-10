library(monocle3)
library(SeuratWrappers)
library(Seurat)

# load data
ExN_sub = readRDS(all,'data/AMY/ExN_subset.rds')

# transfer data
ExN_sub.monocle <- as.cell_data_set(ExN_sub)
ExN_sub.monocle <- estimate_size_factors(ExN_sub.monocle)

# get trajectory
ExN_sub.monocle <- cluster_cells(ExN_sub.monocle,resolution=0.0005)
ExN_sub.monocle <- learn_graph(ExN_sub.monocle)

# set root
ExN_sub.monocle <- order_cells(ExN_sub.monocle, root_pr_nodes='Y_47')

# plot
plot_cells(ExN_sub.monocle,color_cells_by = 'pseudotime',show_trajectory_graph = T,label_cell_groups = F,label_branch_points = F,label_principal_points = F,label_leaves = F,label_roots = F)+NoAxes()
