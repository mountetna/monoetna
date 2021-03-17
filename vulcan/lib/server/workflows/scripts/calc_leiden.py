from archimedes.functions.dataflow import output_path, input_path, input_var, json
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('nn_anndata.h5ad'))

# umap
sc.tl.leiden(scdata, key_added = "leiden_1.0")

##### OUTPUT
with open(output_path('leiden.json'), 'w') as output_file:
    json.dump(scdata.obs['leiden_1.0'].tolist(), output_file)
