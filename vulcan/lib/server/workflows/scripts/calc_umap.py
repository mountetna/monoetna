from archimedes.functions.dataflow import output_path, input_path, input_var
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('pca_anndata.h5ad'))
max_pc = int(input_var('max_pc'))

# Neighborhood and umap
sc.pp.neighbors(scdata, n_neighbors=10, n_pcs= max_pc)
sc.tl.umap(scdata)

##### OUTPUT
scdata.write(output_path('umap_anndata.h5ad'))
