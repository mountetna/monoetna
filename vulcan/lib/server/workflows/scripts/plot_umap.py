from archimedes.functions.dataflow import output_path, input_path, input_var, json, input_json
from archimedes.functions.utils import re
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.plotting import px, pio, colors
from archimedes.functions.magma import connect, question
from archimedes.functions.list import flatten


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
        *value
    ], strip_identifiers=False))

    return [ values.get(id, None) for id in ids ]

pdat = input_json("project_data")[project_name]
color_options = pdat['color_options']

if color_by == 'Cluster':
    color = leiden
elif color_by == 'Tube':
    color = scdata.obs[ 'Record_ID' ]
elif color_by in scdata.raw.var_names:
    color = flatten(scdata.raw.X[ : , scdata.raw.var_names == color_by ].toarray())
elif color_by in color_options.keys():
    color = get(scdata.obs[ 'Record_ID' ], buildTargetPath(color_options[color_by], pdat))
else:
    color = None

##### OUTPUT
fig = px.scatter(
    scdata.obsm['X_umap'], x=0, y=1,
    color_discrete_sequence=colors,
    color=color)

fig.update_layout(
    xaxis_title='UMAP0',
    yaxis_title='UMAP1',
    legend_title=color_by
)
fig.update_coloraxes(colorbar_title_text=color_by)

fig.update_traces(marker={'size': 5, 'opacity': 0.5})

with open(output_path('umap.plotly.json'), 'w') as output_file:
    json.dump(json.loads(pio.to_json(fig)), output_file)
