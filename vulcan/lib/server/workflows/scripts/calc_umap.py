from archimedes.functions.dataflow import output_path, input_path
from archimedes import scanpy

adata = open(input_path('data'), 'r').read()
max_pc = open(input_path('max_pc'), 'r').read()

# Neighborhood and umap
sc.pp.neighbors(adata, n_neighbors=10, n_pcs= max_pc)
sc.tl.umap(adata)

##### OUTPUT
adata.write(output_path('umap_data'))
