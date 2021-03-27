from archimedes.functions.dataflow import output_path, input_path, input_var
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('nn_anndata.h5ad'))
min_dist = float(input_var('min_dist'))
spread = float(input_var('spread'))
num_iters = int(input_var('num_iters'))

# umap
sc.tl.umap(scdata, min_dist = min_dist, spread = spread, maxiter = num_iters)

##### OUTPUT
scdata.write(output_path('umap_anndata.h5ad'))
