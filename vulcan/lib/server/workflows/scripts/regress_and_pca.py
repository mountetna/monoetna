from archimedes.functions.dataflow import output_path, input_path, input_bool
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('normed_anndata.h5ad'))
regress_nCounts = input_bool('regress_counts')
regress_nFeatures = input_bool('regress_genes')
regress_pct_mito = input_bool('regress_pct_mito')
regress_pct_ribo = input_bool('regress_pct_ribo')
regress_tube_id = input_bool('regress_tube_id')

# Regress out and scale
if any([regress_nCounts, regress_nFeatures, regress_pct_mito, regress_pct_ribo, regress_tube_id]):
    regress_on = []
    if regress_nCounts: 
        regress_on.append('n_counts')
    if regress_nFeatures: 
        regress_on.append('n_genes')
    if regress_pct_mito: 
        regress_on.append('pct_counts_mt')
    if regress_pct_ribo: 
        regress_on.append('pct_counts_rb')
    if regress_tube_id: 
        regress_on.append('Record_ID')
    print('No regression is being performed, but the system would otherwise have regressed on:')
    print(regress_on)
    # sc.pp.regress_out(scdata, regress_on)
scdata

sc.pp.scale(scdata, max_value=10)
scdata

# pca
sc.tl.pca(scdata, svd_solver='arpack')

scdata

##### OUTPUT
scdata.write(output_path('pca_anndata.h5ad'))
