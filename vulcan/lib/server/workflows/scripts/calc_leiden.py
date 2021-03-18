from archimedes.functions.dataflow import output_path, input_path, input_var, json
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('nn_anndata.h5ad'))

# umap
leiden_resolution = float(input_var('leiden_resolution'))
sc.tl.leiden(scdata, resolution = leiden_resolution, key_added = "leiden")

##### OUTPUT
with open(output_path('leiden.json'), 'w') as output_file:
    json.dump(scdata.obs['leiden'].tolist(), output_file)
