from archimedes.functions.dataflow import output_path, input_path, input_var
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('nn_anndata.h5ad'))

# umap
sc.tl.umap(scdata)

##### OUTPUT
scdata.write(output_path('umap_anndata.h5ad'))
