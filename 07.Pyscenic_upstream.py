import scanpy as sc
import pandas as pd
import loompy as lp
import numpy as np

# Load data
adata = sc.read("data/AMY/scenic_subset.h5ad")

# Create loom object
row_attrs = {
    "Gene": np.array(Progenitor.var_names) ,
}
col_attrs = {
    "CellID": np.array(Progenitor.obs_names) ,
    "nGene": np.array( np.sum(Progenitor.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(Progenitor.X.transpose() , axis=0)).flatten()
}
lp.create( 'data/AMY/scenic_subset.loom', Progenitor.X.transpose(), row_attrs, col_attrs)

# Build gene regulatory network
f_tfs = "software/SCENIC/hs_hgnc_curated_tfs.txt"
f_loom_path_scenic = "data/AMY/scenic_subset.loom"
!pyscenic grn {f_loom_path_scenic} {f_tfs} -o adj.csv --num_workers 10

# Infer regulon
adjacencies = pd.read_csv("./adj.csv", index_col=False)
f_db_names = 'software/SCENIC/databases/homo_sapiens/hg38/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather'
f_motif_path ='software/SCENIC/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
!pyscenic ctx adj.csv \
    {f_db_names} \
    --annotations_fname {f_motif_path} \
    --expression_mtx_fname {f_loom_path_scenic} \
    --output reg.csv \
    --mask_dropouts \
    --mode "dask_multiprocessing" \
    --num_workers 10
	
# Compute regulon activity score
f_pyscenic_output = "pyscenic_output.loom"
!pyscenic aucell \
    {f_loom_path_scenic} \
    reg.csv \
    --output {f_pyscenic_output} \
    --num_workers 20