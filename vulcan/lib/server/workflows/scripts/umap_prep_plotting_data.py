from archimedes.functions.dataflow import output_path, output_json, input_path, input_var, input_json, buildTargetPath
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.magma import connect, question
from archimedes.functions.list import flatten
from archimedes.functions.environment import project_name
from archimedes.functions.utils import pandas as pd
from archimedes.functions.utils import csv

# Read inputs
scdata = sc.read(input_path('umap_anndata.h5ad'))

leiden = list(map(str, input_json('leiden.json')))

color_by = input_var('color_by')

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

pdat = input_json("project_data")[project_name]
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
plot_presets = {
    'x_by': '0',
    'y_by': '1',
    'color_by': color_by,
    'xlab': 'UMAP_1',
    'ylab': 'UMAP_2'
}

if custom_tooltip:
    dat[hover_name] = hover_text
    plot_presets['hover_data'] = [hover_name]

# Output data_frame and preset options for plotting
dat.to_json(output_path("data_frame"))
output_json(plot_presets, "preset")
