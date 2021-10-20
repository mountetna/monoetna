from archimedes.functions.dataflow import input_path, output_path
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('scdata.h5ad'))
metadata = scdata.obs

metadata.to_json(output_path("metadata"))