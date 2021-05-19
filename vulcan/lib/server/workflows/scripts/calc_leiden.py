from archimedes.functions.dataflow import output_path, input_path, input_var, input_bool, output_json
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('nn_anndata.h5ad'))
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
        ] for clust in range(max(clusts)+1)
    )
    'blank_annots.json'
)

### Additional
# color options
pdat = input_json("project_data")[project_name]

color_options = {
    'Gene': dict([ [ gene_id, None ] for gene_id in scdata.raw.var_names ]),
    'Cluster': None,
    'Manual Annotations': None,
    'Tube': None
}
color_options.update([ [ label, None ] for label in pdat['color_options'].keys() ])

output_json(color_options, 'color_options')