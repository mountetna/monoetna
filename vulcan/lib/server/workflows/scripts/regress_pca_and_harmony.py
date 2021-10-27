from archimedes.functions.dataflow import output_path, input_path, input_bool, input_var
from archimedes.functions.scanpy import scanpy as sc

def output_var(data, name):
    with open(output_path(name), 'w') as output_file:
        output_file.write(data)

scdata = sc.read(input_path('normed_anndata.h5ad'))
regress_nCounts = input_bool('regress_counts')
regress_nFeatures = input_bool('regress_genes')
regress_pct_mito = input_bool('regress_pct_mito')
regress_pct_ribo = input_bool('regress_pct_ribo')
no_batch_string = input_var('no_batch_string')
batch_by = input_var('batch_by')

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
pca_use = 'X_pca'

# harmony batch correction
if not batch_by==no_batch_string:
    sc.external.pp.harmony_integrate(scdata, key=batch_by)
    pca_use='X_pca_harmony'

##### OUTPUT
scdata.write(output_path('pca_anndata.h5ad'))
output_var(pca_use, 'pca_use')
