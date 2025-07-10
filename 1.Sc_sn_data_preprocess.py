# Load packages
import scanpy as sc
import scanpy.external as sce
import os
import matplotlib.pyplot as plt
import time
import re

# Some settings
sc.settings.verbosity=3
sc.set_figure_params(dpi=100,dpi_save=300)
dir_path = '~/AMY/AMY_fetal/CellRanger_pre/'
all_files = [i for i in os.listdir(dir_path) if 'outs_' in i]

# Read data and add metadata
lane_list=['GW08_C','GW09_C','GW09J_N-1','GW09J_N-2','GW10_C','GW11_C-1','GW12_C1','GW12_C2','GW12J_C3','GW13_C','GW13J_N-1','GW13J_N-2','GW13.5J_C-1','GW13.5J_C-2','GW14_C-1','GW14_C-2','GW14_C-3','GW17_N1','GW17_N2','GW19J_N1','GW19J_N2','GW21J_C','GW21J_N-1','GW21J_N-2','GW22J_N1','GW22J_N2','GW22J_N3','GW24_C1','GW24_C2-1','GW24_C2-2','GW24J_N-1','GW24J_N-2','GW25J_C','GW25J_N-1','GW25J_N-2']
week_list=['GW08','GW09','GW09','GW09','GW10','GW11','GW12','GW12','GW12','GW13','GW13','GW13','GW13','GW13','GW14','GW14','GW14','GW17','GW17','GW19','GW19','GW21','GW21','GW21','GW22','GW22','GW22','GW24','GW24','GW24','GW24','GW24','GW25','GW25','GW25']
sample_list=['GW08','GW09#1','GW09J#2','GW09J#2','GW10','GW11','GW12#1','GW12#2','GW12J#3','GW13#1','GW13J#3','GW13J#3','GW13.5J','GW13.5J','GW14','GW14','GW14','GW17','GW17','GW19J','GW19J','GW21J','GW21J','GW21J','GW22J','GW22J','GW22J','GW24#1','GW24#2','GW24#2','GW24J#3','GW24J#3','GW25J','GW25J','GW25J']
sex_list=['female','female','female','female','female','female','female','male','male','male','female','female','male','male','female','female','female','male','male','female','female','male','male','male','male','male','male','female','female','female','male','male','male','male','male']
info_list = ['outs_Amy8GW_20191227_CAB_7000','outs_Amy9GW_20200610_CA_7000','outs_Amy9JGW_20210705_NAB-1_pre_12918','outs_Amy9JGW_20210705_NAB-2_pre','outs_Amy10GW_20191218_CA_7000','outs_Amy11GW_20201218_CA-1_pre','outs_Amy12GW_20200414_CB','outs_Amy12GW_20200603_CA_7000','outs_Amy12JGW_20200813_CA_6056','outs_Amy13GW_20200526_CA_7000','outs_Amy13JGW_20210106_NA-1_pre','outs_Amy13JGW_20210106_NA-2_pre','outs_Amy135JGW_20210105_CA-1_pre','outs_Amy135JGW_20210105_CA-2_pre','outs_Amy14GW_20201209_CA-1_pre','outs_Amy14GW_20201209_CA-2_pre','outs_Amy14GW_20201209_CA-3_pre','outs_Amy17GW_20210331_NA1_pre','outs_Amy17GW_20210331_NA2_pre','outs_Amy19JGW_20210414_NA1_pre','outs_Amy19JGW_20210414_NA2_pre','outs_Amy21JGW_20201109_CA','outs_Amy21JGW_20201109_NA-1_pre','outs_Amy21JGW_20201109_NA-2_pre','outs_Amy22JGW_20210429_NA1_pre','outs_Amy22JGW_20210429_NA2_pre','outs_Amy22JGW_20210429_NA3_pre','outs_Amy24GW_20201203_CA1_11803_pre','outs_Amy24GW_20201203_CA2-1_6517_pre','outs_Amy24GW_20201203_CA2-2_6289_pre','outs_Amy24JGW_20210106_NA-1_pre','outs_Amy24JGW_20210106_NA-2_pre','outs_Amy25JGW_20200809_CA_7750','outs_Amy25JGW_20200809_NA-1_pre','outs_Amy25JGW_20200809_NA-2_pre']
adatas=[]
for i,j in enumerate(all_files):
    temp = sc.read_10x_mtx(dir_path+j+'/outs/filtered_feature_bc_matrix',cache=True,var_names='gene_symbols')
    temp.var_names_make_unique()
    temp.obs['source']=j
    temp.obs['info'] = info_list[i]
    temp.obs_names=[k[:-2] for k in temp.obs_names]
    temp.obs['lane']=lane_list[i]
    temp.obs['week']=week_list[i]
    temp.obs['sample']=sample_list[i]
    temp.obs['sex']=sex_list[i]
    if 'J' in sample_list[i]:
        temp.obs['percise']='Yes'
    else:
        temp.obs['percise']='No'
    adatas.append(temp)
    
# Remove doublets
for i,temp in enumerate(adatas):
    sce.pp.scrublet(temp,expected_doublet_rate=doublet_list[i])   #doublet ratio predicted by no-liner regress
    fig = sce.pl.scrublet_score_distribution(temp,return_fig=True)
    fig.savefig('./1.scrublet/doublet_'+all_files[i]+'.png')
    
# Merge objects
raw=adatas[0].concatenate(adatas[1:],batch_key=None,batch_categories=[str(i+1) for i in range(len(adatas))])
raw.write('data/AMY/raw.h5ad')

# QC
adata=raw[~raw.obs['predicted_doublet'],:]
sc.pp.filter_genes(adata,min_cells=27)  # cells / 10000
adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
adata.var['hb']=adata.var_names.str.contains(("^HB[ABDEGMQZ]\d?$"))
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','hb','ribo'], percent_top=None, log1p=False, inplace=True)
sc.pl.highest_expr_genes(adata, n_top=20)

# Filter cells
print(adata.obs.apply(lambda x:(x['n_genes_by_counts']<8000) & (x['n_genes_by_counts']>800),axis=1).value_counts())
adata=adata[adata.obs.apply(lambda x:(x['n_genes_by_counts']<8000) & (x['n_genes_by_counts']>800),axis=1),:]
print(adata.obs.apply(lambda x:x['pct_counts_mt']<5,axis=1).value_counts())
adata=adata[adata.obs.apply(lambda x:x['pct_counts_mt']<5,axis=1),:]
print(adata.obs.apply(lambda x:x['pct_counts_hb']<0.2,axis=1).value_counts())
adata=adata[adata.obs.apply(lambda x:x['pct_counts_hb']<0.2,axis=1),:]
print(adata.obs.apply(lambda x:x['pct_counts_ribo']<20,axis=1).value_counts())
adata=adata[adata.obs.apply(lambda x:x['pct_counts_ribo']<20,axis=1),:]
print(adata.obs.apply(lambda x:x['total_counts']<60000,axis=1).value_counts())
adata=adata[adata.obs.apply(lambda x:x['total_counts']<60000,axis=1),:]

# Visualization
fig = plt.figure(figsize=(40,36))
grid = plt.GridSpec(nrows=5, ncols=4, figure=fig)
ax1=plt.subplot(grid[0,1:4])
sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4, multi_panel=False,groupby=('lane'), rotation=90,size=0.5,ax=ax1,show=False,order=lane_list)
ax2=plt.subplot(grid[0,0])
sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4, multi_panel=False,ax=ax2,show=False,size=0.8)
ax3=plt.subplot(grid[1,1:4])
sc.pl.violin(adata, ['total_counts'], jitter=0.4, multi_panel=False,groupby=('lane'), rotation=90,size=0.5,ax=ax3,show=False,order=lane_list)
ax4=plt.subplot(grid[1,0])
sc.pl.violin(adata, ['total_counts'], jitter=0.4, multi_panel=False,ax=ax4,show=False,size=0.8)
ax5=plt.subplot(grid[2,1:4])
sc.pl.violin(adata, ['pct_counts_mt'], jitter=0.4, multi_panel=False,groupby=('lane'), rotation=90,size=0.5,ax=ax5,show=False,order=lane_list)
ax6=plt.subplot(grid[2,0])
sc.pl.violin(adata, ['pct_counts_mt'], jitter=0.4, multi_panel=False,ax=ax6,show=False,size=0.8)
ax7=plt.subplot(grid[3,1:4])
sc.pl.violin(adata, ['pct_counts_hb'], jitter=0.4, multi_panel=False,groupby=('lane'), rotation=90,size=0.5,ax=ax7,show=False,order=lane_list)
ax8=plt.subplot(grid[3,0])
sc.pl.violin(adata, ['pct_counts_hb'], jitter=0.4, multi_panel=False,ax=ax8,show=False,size=0.8)
ax9=plt.subplot(grid[4,1:4])
sc.pl.violin(adata, ['pct_counts_ribo'], jitter=0.4, multi_panel=False,groupby=('lane'), rotation=90,size=0.5,ax=ax9,show=False,order=lane_list)
ax10=plt.subplot(grid[4,0])
sc.pl.violin(adata, ['pct_counts_ribo'], jitter=0.4, multi_panel=False,ax=ax10,show=False,size=0.8)
fig.subplots_adjust(hspace=0.5)
fig.savefig('./2.QC/QC_after.png',dpi=300,bbox_inches='tight')

# Normalized data
adata.layers['counts']=adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers['log']=adata.X.copy()
adata.raw=adata
sc.pp.highly_variable_genes(adata)
adata = adata[:,adata.var.highly_variable]
sc.pp.scale(adata)

# Save data
adata.write('data/AMY/afterQC.h5ad')