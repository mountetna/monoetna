from archimedes.functions.dataflow import output_path, output_json, input_path, input_var, input_json, buildTargetPath
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.list import flatten
from archimedes.functions.utils import pandas as pd
from archimedes.functions.utils import re

# Read inputs
scdata = sc.read(input_path('scdata.h5ad'))

color_by = input_var('color_by')

pdat = input_json("project_data")
color_options = pdat['color_options']

# Obtain color data
custom_tooltip = False
if color_by in ['leiden', 'Clustering (leiden)']:
    color = scdata.obs[ 'leiden' ].values
    custom_tooltip = True
    sets = input_json('top10.json')
    texts = dict([
        [
            str(clust),
            ' ' + ', '.join(sets[clust])
        ] for clust in list(sets.keys()) ])
    hover_name = 'top10 markers'
    hover_text = [texts[str(val)] for val in color]
elif color_by in list(scdata.obs.columns):
    color = scdata.obs[ color_by ].values
elif color_by in scdata.raw.var_names:
    color = flatten(scdata.raw.X[ : , scdata.raw.var_names == color_by ].toarray())
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
