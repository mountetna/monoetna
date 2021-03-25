from archimedes.functions.dataflow import output_path, input_path, input_var, json, input_json
from archimedes.functions.utils import re
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.plotting import px, pio, colors
from archimedes.functions.magma import connect, question

scdata = sc.read(input_path('umap_anndata.h5ad'))

leiden = list(map(str, input_json('leiden.json')))

color_by = input_var('color_by')

magma = connect()

def get(ids, value):
    ids = ids if type(ids) == list else ids.tolist()
    values = dict(question(magma, [
        'sc_seq',
        [ '::identifier', '::in', ids ],
        '::all',
        value
    ], strip_identifiers=False))

    return [ values.get(id, None) for id in ids ]

if color_by == 'Cluster':
    color = leiden
elif color_by == 'Experiment':
    color = [
        re.sub( 'PT.[0-9]+$', '', v.split('-')[2])
        for v in scdata.obs[ 'Record_ID']
    ]
elif color_by == 'Tissue':
    color = get(scdata.obs[ 'Record_ID' ], [ 'biospecimen_group', 'biospecimen_type' ])
elif color_by == 'Pool':
    color = get(scdata.obs[ 'Record_ID' ], [ 'sc_seq_pool', '::identifier' ])
elif color_by == 'Biospecimen Group':
    color = get(scdata.obs[ 'Record_ID' ], [ 'biospecimen_group', '::identifier' ])
elif color_by in scdata.var_names:
    color = scdata.X[ : , scdata.var_names == color_by ]
else:
    color = None

##### OUTPUT
fig = px.scatter(
    scdata.obsm['X_umap'], x=0, y=1,
    color_discrete_sequence=colors,
    color=color)

with open(output_path('umap.plotly.json'), 'w') as output_file:
    json.dump(json.loads(pio.to_json(fig)), output_file)
