from scipy import sparse
import scanpy as sc
import anndata as ad
import pandas as pd
import os
os.environ["R_HOME"] = "./Software/miniconda3/envs/R_main/lib/R"
import scanpy as sc
import STAGATE as sta
import tensorflow as tf
import sys
tf.compat.v1.disable_eager_execution()

# set datadir
gem_path = "./DATA/sterseq/Human Spatial/"
out_dir = './AMY_spatial/result/'

# read gem file
for index,i in enumerate(os.listdir(gem_path)):
	df = pd.read_csv(
    i,
    sep="\t",
    comment="#"
	)

	# substract infomation
	expr = (
		df.groupby(["CellID", "geneID"])["MIDCount"]
		.sum()
		.unstack(fill_value=0)
	)
	coords = (
		df.groupby("CellID")[["x", "y"]]
		.mean()
		.loc[expr.index])

	# create anndata file
	adata = ad.AnnData(expr)
	
	# spatial coords
	adata.obsm["spatial"] = coords.values
	
	# set index
	adata.obs_names = expr.index.astype(str)
	adata.var_names = expr.columns.astype(str)
	adata.X = sparse.csr_matrix(adata.X).copy()

    # normalize data
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    
    # STAGATE pipeline
    sc.pp.highly_variable_genes(adata,n_top_genes=2000)
    sta.Cal_Spatial_Net(adata, k_cutoff=30,model='KNN')
    sta.Stats_Spatial_Net(adata)
    adata = sta.train_STAGATE(adata)
    adata = sta.mclust_R(adata, used_obsm='STAGATE', num_cluster=30)
    adata.write(out_dir+'data_'+str(index)+'h5ad')

