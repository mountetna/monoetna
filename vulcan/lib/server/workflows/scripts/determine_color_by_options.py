from archimedes.functions.dataflow import input_path, output_json, input_json
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.environment import project_name

scdata = sc.read(input_path('leiden_anndata.h5ad'))
pdat = input_json("project_data")[project_name]

color_options = {
    'Gene': dict([ [ gene_id, None ] for gene_id in scdata.raw.var_names ]),
    'Cluster': None,
    'Manual Annotations': None,
    'Tube': None
}
color_options.update([ [ label, None ] for label in pdat['color_options'].keys() ])

output_json(color_options, 'color_options')
