from archimedes.functions.dataflow import output_path, input_path, input_bool
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('normed_anndata.h5ad'))
regress_nCounts = input_bool('regress_counts')
regress_nFeatures = input_bool('regress_genes')
regress_pct_mito = input_bool('regress_pct_mito')
regress_pct_ribo = input_bool('regress_pct_ribo')

# Regress out and scale
if any([regress_nCounts, regress_nFeatures, regress_pct_mito, regress_pct_ribo]):
    regress_on = []
    if regress_nCounts: 
        regress_on.append('n_counts')
    if regress_nFeatures: 
        regress_on.append('n_genes')
    if regress_pct_mito: 
        regress_on.append('pct_counts_mt')
    if regress_pct_ribo: 
        regress_on.append('pct_counts_rb')
    sc.pp.regress_out(scdata, regress_on)

sc.pp.scale(scdata, max_value=10)

# pca
sc.tl.pca(scdata, svd_solver='arpack')

##### OUTPUT
scdata.write(output_path('pca_anndata.h5ad'))
