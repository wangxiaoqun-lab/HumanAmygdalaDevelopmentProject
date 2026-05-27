import scanpy as sc
import scvelo as scv
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt

# load data
adata = sc.read('../DATA/all.h5ad')

# load spliced/unspliced matrix
scv_list = ['GW09_1','GW09_1_re','GW09_2','GW09_3','GW10_1','GW10_2','GW11','GW13','GW13_re','GW13.6_1','GW13.6_2','GW14','GW17','GW18','GW19_part1','GW19_part2','GW21','GW21_re','GW22.5_part1','GW22.5_part2','GW22.5_part3','GW24','GW24_re','GW25_1','GW25_1_re','GW25_2_part1','GW25_2_part2','GW25_2_part3']
for lane in []:
    temp = adata[adata.obs['label']==lane].copy()
    loom = sc.read('../DATA/loomfile/human_amy_sn_loom/'+lane+'.loom')
    loom = scv.utils.merge(temp, loom)
    scv_list.append(loom)

# merge data
adata_velo = scv_list[0].concatenate(scv_list[1:])

# nomalize data
scv.pp.filter_and_normalize(
    adata_velo, min_shared_counts=20, n_top_genes=None, subset_highly_variable=False
)

# get moments
scv.pp.moments(adata_velo, n_pcs=30)

# get velocity graph
scv.tl.velocity_graph(adata_velo)

# plot and save data
scv.pl.velocity_embedding_stream(adata_velo, basis='umap',color='celltype_new')
adata_velo.write('./DATA/adata_velocity_version2.h5ad')