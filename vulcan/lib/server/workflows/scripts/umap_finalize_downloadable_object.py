from archimedes.functions.dataflow import input_path, input_json, output_path
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('scdata.h5ad'))

# Add Manual Annots
annots = input_json('annots.json')
leiden = list(map(str, input_json('leiden.json')))
data = [annots[x] for x in leiden]
scdata.obs[ 'Manual_Annotations' ] = data

# Note: Clusters & Tube are already in the scdata.obs as 'leiden' and 'Record_ID'.
#   So are color-options of the project-definitions!

# Output the anndata object
scdata.write(output_path('umap_workflow_anndata.h5ad'))
