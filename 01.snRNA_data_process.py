# Load packages
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
import time
import re
from sklearn.linear_model import LinearRegression

# Some settings
sc.settings.verbosity=3
sc.set_figure_params(dpi=100,dpi_save=300)
dir_path = '~/AMY/AMY_fetal/'

# Read data
adata = sc.read(dir_path + 'DATA/AMY_all.h5ad')

# QC
## mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
## hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt",  "hb"],percent_top=None,log1p=False,inplace=True)
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb"],
    jitter=0.4,
    multi_panel=True,
)

# Filter cells
group_key = 'label'
counts_key = 'total_counts'
mt_key = 'pct_counts_mt'
hb_key = 'pct_counts_hb'

low_q = 0.01
high_q = 0.98
mt_thresh = 5.0
hb_thresh = 1.0

adata.obs['pass_basic_qc'] = True
for g, idx in adata.obs.groupby(group_key).groups.items():
    counts = adata.obs.loc[idx, counts_key]
    low_cut = counts.quantile(low_q)
    high_cut = counts.quantile(high_q)
    mask = (counts >= low_cut) & (counts <= high_cut)
    adata.obs.loc[idx, 'pass_basic_qc'] &= mask
adata.obs['pass_basic_qc'] &= (
    (adata.obs[mt_key] <= mt_thresh) &
    (adata.obs[hb_key] <= hb_thresh)
)
adata_qc = adata[adata.obs['pass_basic_qc']].copy()

# Visualization
fig,ax = plt.subplots(figsize=(10,5))
sc.pl.violin(adata_qc, ['total_counts'], groupby='label',multi_panel=False,stripplot=False,rotation=90,ax=ax)

# Remove doublet
sce.pp.scrublet(adata_qc,batch_key='label')
sce.pl.scrublet_score_distribution(adata_qc)
adata_qc.obs['predicted_doublet'] = adata_qc.obs['doublet_score'] > 0.2

# Check R2
adata_final = adata_qc[~ adata_qc.obs['predicted_doublet']].copy()
batches = adata_final.obs['label'].unique()
adata_final.obs['R2_outlier'] = False
percentile = 0.01
for i, batch in enumerate(batches):
    subset = adata_final.obs[adata_final.obs['label'] == batch]
    X = subset['total_counts'].values.reshape(-1,1)
    y = subset['n_genes_by_counts'].values
    model = LinearRegression().fit(X, y)
    y_pred = model.predict(X)
    r2 = model.score(X, y)
    residuals = y - y_pred
    lower_thresh = np.percentile(residuals, percentile*100)
    upper_thresh = np.percentile(residuals, 100 - percentile*100)
    is_outlier = (residuals < lower_thresh) | (residuals > upper_thresh)
    adata_final.obs.loc[subset.index, 'R2_outlier'] = is_outlier
adata = adata_final[~ adata_final.obs['R2_outlier']]

# Save data
adata.write(dir_path +'DATA/AMY_all_after_QC.h5ad')

# normalize
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)

# dim reduction
sc.pp.highly_variable_genes(adata,n_top_genes=2000,batch_key='label')
sc.pp.pca(adata)
sc.pl.pca_variance_ratio(adata,n_pcs=50,log=True)
sce.pp.bbknn(adata,neighbors_within_batch=10,batch_key='Sample',n_pcs=30)
sc.tl.umap(adata,min_dist=0.8)
sc.pl.umap(adata,color='label')


# clustering and annotation
sc.tl.leiden(adata,resolution=1,key_added='leiden_new')
## re-cluster MIX cells
sc.tl.leiden(adata,resolution=0.1,restrict_to=('leiden_new',['22']),key_added='leiden_R_new')
adata.obs['celltype_new'] = adata.obs['leiden_R_new'].replace(['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22,0','22,1'],
                                                              ['InN','InN','Progenitor','InN','ExN','InN','ExN','InN','ExN','ExN','ExN','ExN','InN','InN','Astrocyte','ExN','InN','ExN','ExN','ExN','OPC&OL','InN','Microglia','Endo&Peri'])

# marker genes
fig = sc.pl.umap(adata,color=['NES','VIM','TOP2A','GAD1','GAD2','SLC17A6','NEUROD2','GFAP','AQP4','PDGFRA','OLIG2','PTPRC','C3','CLDN5','FLT1'],vmax='p99.9',size=1.5,cmap='Reds',ncols=5,return_fig=True)
fig.savefig('./Figs/UMAP_marker_genes.pdf',dpi=400,bbox_inches='tight')

#compute DEG
sc.tl.rank_genes_groups(adata,groupby='celltype_new',n_genes=200,method='wilcoxon')

# save data
adata.write(dir_path + 'DATA/all.h5ad')