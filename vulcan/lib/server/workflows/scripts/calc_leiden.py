from archimedes.functions.dataflow import output_path, input_path, input_var, input_bool, output_json
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('umap_anndata.h5ad'))
use_weights = input_bool('use_weights')

# Calculate leiden clustering
leiden_resolution = float(input_var('leiden_resolution'))
sc.tl.leiden(scdata, resolution = leiden_resolution, key_added = "leiden", use_weights = use_weights)

### Output
scdata.write(output_path('leiden_anndata.h5ad'))

clusts = scdata.obs['leiden'].tolist()

output_json(clusts, 'leiden.json')
output_json(
    dict(
        [
            str(clust),
            str(clust)
        ] for clust in range(max(
            [int(str) for str in clusts]
            )+1)
    ),
    'blank_annots.json'
)