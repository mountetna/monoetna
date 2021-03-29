from archimedes.functions.dataflow import output_path, output_json, input_path, input_var
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('nn_anndata.h5ad'))
min_dist = float(input_var('min_dist'))
spread = float(input_var('spread'))
num_iters = int(input_var('num_iters'))

# umap
sc.tl.umap(scdata, min_dist = min_dist, spread = spread, maxiter = num_iters)

# Bring pre-subsetting/regression/scaling expression data back to forefront
scdata.raw.to_adata()

# color options
color_options = {
    'Cluster': None,
    'Experiment': None,
    'Tissue': None,
    'Pool': None,
    'Biospecimen Group': None,
    'Gene': dict([ [ gene_id, None ] for gene_id in scdata.var_names ])
}

##### OUTPUT
scdata.write(output_path('umap_anndata.h5ad'))

output_json(color_options, 'color_options')

