from archimedes.functions.dataflow import input_path, output_json
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('scdata.h5ad'))
metadata = scdata.obs

keys = metadata.columns

discrete_metas = list()
continuous_metas = list()

for key in keys:
    if all(map(lambda x: isinstance(x, (str,bool)), metadata[key])):
        discrete_metas.append(key)
    else:
        continuous_metas.append(key)

def blank_dict(keys):
    return dict([k, None] for k in keys)

output_json(discrete_metas, "discrete_metas")
output_json(continuous_metas, "continuous_metas")