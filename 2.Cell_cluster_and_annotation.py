import scanpy as sc
import scanpy.external as sce
import os
import matplotlib.pyplot as plt
import time
import re

# load data
adata = sc.read('data/AMY/afterQC.h5ad')

#maker genes
classical_markers={'All Neurons':['SNAP25','MAP2', 'SYT1', 'SYP', 'RBFOX3','KIF5C'],
                   'Excitatory Neuron':['SLC17A7', 'NRGN','SLC17A6','CAMK2A','NEUROD2'],
                   'Inhibitory Neuron|Interneuron':['GAD1', 'GAD2', 'SLC32A1', 'SLC6A1'],
                   'Astrocyte':['AQP4', 'GFAP', 'FGFR3', 'SLC1A2', 'SLC1A3', 'ALDH1A1', 'ALDH1L1','SLC7A10','LFNG'],
                   'Oligodendrocyte':['MOG', 'MOBP', 'MBP', 'PLP1', 'OPALIN','HAPLN2'],
                   'Oligodendrocyte Precursor cell':['PDGFRA', 'CSPG4', 'VCAN', 'BCAN', 'PCDH15', 'APOD','GPR17', 'OLIG2'],
                   'Microglia':['PTPRC', 'P2RY12', 'C3', 'HLA-DRA', 'CD74', 'CX3CR1', 'AIF1','C1QC','TMEM119'],
                   'Endothelial cell':['PECAM1', 'CLDN5', 'VWF', 'NOSTRIN', 'FLT1'],
                   'Pericyte':['DCN', 'COL1A2', 'PDGFRB', 'TAGLN','VTN','ENKUR','CCDC153','TBX18','MYH11','RGS5'],
                   'Progenitor':['ASCL1', 'MKI67','SOX2','PAX6','NES','VIM','EOMES','HOPX','NKX2-1']
                  }
cell_cycle_genes=['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7','DTL','PRIM1','UHRF1','MLF1IP','HELLS','RFC2','RPA2','NASP','RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2','RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2','USP1','CLSPN','POLA1','CHAF1B','BRIP1','E2F8','HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2','NUF2','CKS1B','MKI67','TMPO','CENPF','TACC3','FAM64A','SMC4','CCNB2','CKAP2L','CKAP2','AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1','KIF20B','HJURP','CDCA3','HN1','CDC20','TTK','CDC25C','KIF2C','RANGAP1','NCAPD2','DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1','ANLN','LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA']
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

#compute cell cycle score
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# remove batch effert
sc.tl.pca(adata,n_comps=100)
sce.pp.bbknn(adata,batch_key='sample',n_pcs=40)

#dim reduction
sc.tl.umap(adata)

#cluster cells
for i in [1,1.5,2,2.5]:
    sc.tl.leiden(adata,resolution=i,key_added='res_'+str(i))
for i in ['res_1','res_1.5','res_2','res_2.5']:
    sc.pl.umap(adata,color=i)
    
# annotation cells
cluster_num=adata.obs['res_1'].dtype.categories.tolist()
cluster_type=['ExN','InN','ExN','Progenitor','ExN','Progenitor','MIX','ExN','InN','ExN','ExN','MIX','InN','InN','InN','ExN','ExN','Astrocyte','InN','MIX','InN','ExN','ExN','InN','OPC&OL','Microglia','InN','Endo&Peri']
adata.obs['cell_type']=adata.obs['res_1'].replace(cluster_num,cluster_type)
sc.pl.umap(adata, color=['cell_type'],use_raw=False)
sc.pl.umap(adata, color=['phase'],use_raw=False)
sc.pl.umap(adata, color=['percise'],use_raw=False)

# re-cluster MIX cells
sc.tl.leiden(adata,resolution=2,restrict_to=('cell_type',['MIX']),key_added='temp')
cluster_num=[i for i in adata.obs['temp'].cat.categories if 'MIX' in i]
cluster_type=['InN','ExN','ExN','MIX','InN','MIX','MIX','InN','InN','InN','ExN','MIX','ExN','InN','ExN','InN','ExN','ExN','InN','InN','InN','ExN','InN','ExN','MIX','MIX']
adata.obs['cell_type']=adata.obs['temp'].replace(cluster_num,cluster_type)
sc.pl.umap(adata, color=['cell_type'],use_raw=False)
sc.pl.umap(adata, color=['phase'],use_raw=False)
sc.pl.umap(adata, color=['percise'],use_raw=False)

# third round cluster
sc.tl.leiden(adata,resolution=2,restrict_to=('cell_type',['MIX']),key_added='temp')
cluster_num=[i for i in adata.obs['temp'].cat.categories if 'MIX' in i]
cluster_type=['ExN','ExN','ExN','ExN','ExN','InN','ExN','ExN','ExN','ExN','InN','ExN','ExN','ExN','InN','ExN']
adata.obs['cell_type']=adata.obs['temp'].replace(cluster_num,cluster_type)
sc.pl.umap(adata, color=['cell_type'],use_raw=False)
sc.pl.umap(adata, color=['phase'],use_raw=False)
sc.pl.umap(adata, color=['percise'],use_raw=False)

# re-index
list1 = adata.obs['leiden_R'].cat.categories
list2 = ['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50']
adata.obs['leiden'] = adata.obs['leiden_R'].replace(list1,list2)

#compute DEG
sc.tl.rank_genes_groups(adata,groupby='cell_type',n_genes=200,method='wilcoxon')
# save data
adata.write('data/AMY/result.h5ad')