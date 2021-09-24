from archimedes.functions.dataflow import input_path, input_json, input_var, output_path
from archimedes.functions.scanpy import scanpy as sc

# Load
scdata = sc.read(input_path('scdata.h5ad'))
vals = list(input_json('selected_vals'))
meta = input_var('meta')

# Subset & Output
scdata[scdata.obs[meta].isin(vals)].write(output_path('scdata.h5ad'))