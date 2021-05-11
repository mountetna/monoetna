from archimedes.functions.dataflow import output_path, input_path, input_var
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('pca_anndata.h5ad'))

# Neighborhood
max_pc = int(input_var('max_pc'))
n_neighbors = int(input_var('n_neighbors'))
sc.pp.neighbors(scdata, n_neighbors=n_neighbors, n_pcs= max_pc)

##### OUTPUT
scdata.write(output_path('nn_anndata.h5ad'))
