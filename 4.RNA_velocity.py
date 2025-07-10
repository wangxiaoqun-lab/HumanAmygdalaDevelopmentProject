import scanpy as sc
import scvelo as scv
import scanpy.external as sce
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

scv.settings.verbosity=3
sc.settings.verbosity = 3

# loom files get from velocyto
loom_dir = 'data/loom_files/output_looms/'
source = os.listdir(loom_dir)

# load data
ldatas=[]
rdata = sc.read('data/AMY/afterQC.h5ad')
for i,j in enumerate(source):
    temp = scv.read_loom(loom_dir+j+'.loom')
    temp.obs_names = [k.split(':')[-1][0:-1]+'-'+str(i+1) for k in temp.obs_names]
    temp.var_names_make_unique()
    ldatas.append(temp)
ldata = ldatas[0].concatenate(ldatas[1:],batch_key=None,batch_categories=None,index_unique=None)

# merge data
adata = scv.utils.merge(rdata,ldata[rdata.obs_names,rdata.var_names])

# nolmalize
scv.pp.normalize_per_cell(adata)
scv.pp.log1p(adata)

# compute RNA velocity
scv.pp.moments(adata)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata,n_jobs=6)

# plot
scv.pl.velocity_embedding_grid(adata,density=0.8,arrow_size=(20,10,20),color='cell_type',arrow_length=6,smooth=True,groups=['Progenitor','InN','ExN','Astrocyte','OPC&OL'],min_mass=10,arrow_color='black')