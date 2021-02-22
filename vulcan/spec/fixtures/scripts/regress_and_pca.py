from archimedes.functions.dataflow import output_path, input_path
from archimedes import scanpy

adata = open(input_path('data'), 'r').read()
regress_nCounts = open(input_path('counts'), 'r').read()
regress_nFeatures = open(input_path('genes'), 'r').read()
regress_per_mito = open(input_path('per_mito'), 'r').read()
regress_per_ribo = open(input_path('per_ribo'), 'r').read()

# Regress out and scale
regress_on = ['total_counts', 'total_genes', 'pct_counts_mt', 'pct_counts_rb']
regress_on = regress_on[regress_nCounts, regress_nFeatures, regress_per_mito, regress_per_ribo]
sc.pp.regress_out(adata, regress_on)
sc.pp.scale(adata, max_value=10)

# pca
sc.tl.pca(adata, svd_solver='arpack')

##### OUTPUT
adata.write(output_path('pca_output_data'))
