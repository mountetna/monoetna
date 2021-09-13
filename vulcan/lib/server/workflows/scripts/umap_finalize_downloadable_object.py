from archimedes.functions.dataflow import input_path, output_json, input_json, output_path, buildTargetPath
from archimedes.functions.magma import connect, question
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.environment import project_name
from archimedes.functions.utils import re

scdata = sc.read(input_path('scdata.h5ad'))

# Prep magma querying
magma = connect()
def get(ids, value):
    ids = ids if type(ids) == list else ids.tolist()
    values = dict(question(magma, [
        'sc_seq',
        [ '::identifier', '::in', ids ],
        '::all',
        *value
    ], strip_identifiers=False))

    return [ values.get(id, None) for id in ids ]

# Add data for all color_by options
pdat = input_json("project_data")[project_name]
color_options = pdat['color_options']

for color_by in color_options.keys():
    data = get(scdata.obs[ 'Record_ID' ].values, buildTargetPath(color_options[color_by], pdat))
    label = re.sub(" ", "_", color_by)
    scdata.obs[ label ] = data

# Add Manual Annots
annots = input_json('annots.json')
leiden = list(map(str, input_json('leiden.json')))
data = [annots[x] for x in leiden]
scdata.obs[ 'Manual_Annotations' ] = data

# Note: Clusters & Tube are already in the scdata.obs as 'leiden' and 'Record_ID'.

# Output the anndata object
scdata.write(output_path('umap_workflow_anndata.h5ad'))
