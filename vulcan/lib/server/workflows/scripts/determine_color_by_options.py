from archimedes.functions.dataflow import input_path, output_json
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('leiden_anndata.h5ad'))

color_options = {
    'Gene': dict([ [ gene_id, None ] for gene_id in scdata.raw.var_names ]),
    'Clustering (leiden)': None,
    # Because this is not added into to the leiden_anndata.h5ad
    'Manual_Annotations': None
}
# All color-options will have already been compiled into the object's metadata!
color_options.update([ [ label, None ] for label in list(scdata.obs.columns) ])

output_json(color_options, 'color_options')
