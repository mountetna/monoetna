from archimedes.functions.dataflow import output_path, input_path, input_var, json, input_json, buildTargetPath
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.plotting import pio, scatter_plotly
from archimedes.functions.magma import connect, question
from archimedes.functions.list import flatten
from archimedes.functions.utils import pandas as pd

# Read inputs
scdata = sc.read(input_path('umap_anndata.h5ad'))

leiden = list(map(str, input_json('leiden.json')))

color_by = input_var('color_by')

px_args = {}

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

pdat = input_json("project_data")
color_options = pdat['color_options']

# Obtain color data
custom_tooltip = False
if color_by == 'Cluster':
    color = leiden
    custom_tooltip = True
    sets = input_json('top10.json')
    texts = dict([
        [
            str(clust),
            ' ' + ', '.join(sets[clust])
        ] for clust in list(sets.keys()) ])
    hover_name = 'top10 markers'
    hover_text = [texts[str(val)] for val in leiden]
elif color_by == 'Manual Annotations':
    annots = input_json('annots.json')
    color = [annots[x] for x in leiden]
elif color_by == 'Tube':
    color = scdata.obs[ 'Record_ID' ].values
elif color_by in scdata.raw.var_names:
    color = flatten(scdata.raw.X[ : , scdata.raw.var_names == color_by ].toarray())
elif color_by in color_options.keys():
    color = get(scdata.obs[ 'Record_ID' ].values, buildTargetPath(color_options[color_by], pdat))
else:
    color = None

# Make data.frame and set extra plot settings
dat = pd.DataFrame(scdata.obsm['X_umap'])
dat[color_by] = color

if custom_tooltip:
    dat[hover_name] = hover_text
    px_args['hover_data'] = [hover_name]

# Make & output plot
fig = scatter_plotly(
    dat, 0, 1, color_by,
    px_args = px_args,
    xlab = 'UMAP1',
    ylab = 'UMAP2',
    color_order = 'increasing',
    order_when_continuous_color = True,
    size = 5)

with open(output_path('umap.plotly.json'), 'w') as output_file:
    json.dump(json.loads(pio.to_json(fig)), output_file)
