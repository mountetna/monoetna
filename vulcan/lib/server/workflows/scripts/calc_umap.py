from archimedes.functions.dataflow import output_path, output_json, input_path, input_var
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('nn_anndata.h5ad'))
min_dist = float(input_var('min_dist'))
spread = float(input_var('spread'))
num_iters = int(input_var('num_iters'))

# umap
sc.tl.umap(scdata, min_dist = min_dist, spread = spread, maxiter = num_iters)

# color options
pdat = input_json("project_data")[project_name]

color_options = {
    'Gene': dict([ [ gene_id, None ] for gene_id in scdata.raw.var_names ]),
    'Cluster': None,
    'Tube': None
}
color_options.update([ [ label, None ] for label in pdat['color_options'].keys() ])

##### OUTPUT
scdata.write(output_path('umap_anndata.h5ad'))

output_json(color_options, 'color_options')

