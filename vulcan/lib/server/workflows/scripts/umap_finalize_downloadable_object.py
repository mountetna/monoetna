from archimedes.functions.dataflow import input_path, output_json, input_json, output_path, buildTargetPath, parseModelAttr
from archimedes.functions.magma import connect, question
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.environment import project_name
from archimedes.functions.utils import re
from archimedes.functions.list import unique

scdata = sc.read(input_path('scdata.h5ad'))
pdat = input_json("project_data")[project_name]

# Prep magma querying
seq_model = parseModelAttr(pdat['seq_h5_counts_data'])['model']
magma = connect()
def get(ids, value):
    ids = ids if type(ids) == list else ids.tolist()
    values = dict(question(magma, [
        seq_model,
        [ '::identifier', '::in', ids ],
        '::all',
        *value
    ], strip_identifiers=False))

    return [ values.get(id, None) for id in ids ]

# Add data for all color_by options
color_options = pdat['color_options']

recs = scdata.obs[ 'Record_ID' ].values
for color_by in color_options.keys():
    rec_data = dict(map(
        lambda x: [x, get([x], buildTargetPath(color_options[color_by], pdat))[0]],
        unique(recs) ))
    label = re.sub(" ", "_", color_by)
    #raise Exception(rec_data)
    scdata.obs[ label ] = list( map( lambda x: rec_data[x], recs) )

# Add Manual Annots
annots = input_json('annots.json')
leiden = list(map(str, input_json('leiden.json')))
data = [annots[x] for x in leiden]
scdata.obs[ 'Manual_Annotations' ] = data

# Note: Clusters & Tube are already in the scdata.obs as 'leiden' and 'Record_ID'.

# Output the anndata object
scdata.write(output_path('umap_workflow_anndata.h5ad'))
