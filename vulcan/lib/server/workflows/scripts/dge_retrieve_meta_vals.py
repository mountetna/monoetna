from archimedes.functions.dataflow import input_path, output_json, input_var
from archimedes.functions.list import unique
from archimedes.functions.scanpy import scanpy as sc

metadata = sc.read(input_path('scdata.h5ad')).obs
meta = input_var('meta')
output_json(unique(list(metadata[meta])), 'opts')
