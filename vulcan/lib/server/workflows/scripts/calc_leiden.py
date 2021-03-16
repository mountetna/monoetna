from archimedes.functions.dataflow import output_path, input_path, input_var
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('nn_anndata.h5ad'))

# umap
sc.tl.leiden(scdata, key_added = "leiden_1.0")

##### OUTPUT
scdata.write(output_path('leiden_anndata.h5ad'))
