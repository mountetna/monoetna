from archimedes.functions.dataflow import output_path, input_path, input_var
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('normed_anndata.h5ad'))
regress_nCounts = bool(input_var('regress_counts'))
regress_nFeatures = bool(input_var('regress_genes'))
regress_pct_mito = bool(input_var('regress_pct_mito'))
regress_pct_ribo = bool(input_var('regress_pct_ribo'))

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
    print(regress_on)
    sc.pp.regress_out(scdata, regress_on)
scdata

sc.pp.scale(scdata, max_value=10)
scdata

# pca
sc.tl.pca(scdata, svd_solver='arpack')

scdata

##### OUTPUT
scdata.write(output_path('pca_anndata'))
