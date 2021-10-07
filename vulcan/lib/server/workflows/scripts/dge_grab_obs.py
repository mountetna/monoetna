from archimedes.functions.dataflow import input_path, output_path
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.environment import project_name

data_target = "/app/built" + input_path('data_url').split(project_name)[1]

# scdata = sc.read(input_path('data_url'))
scdata = sc.read(data_target)
metadata = scdata.obs

metadata.to_json(output_path("metadata"))