from archimedes.functions.dataflow import input_json, input_var, input_path, output_json
from archimedes.functions.scanpy import scanpy as sc

setup = input_json('setup')
scdata = sc.read(input_path('scdata.h5ad'))
test_method = input_var('test_method')

output_json(setup, 'diffexp.csv')

# out: [diffexp.csv, filtered_diffexp.csv]