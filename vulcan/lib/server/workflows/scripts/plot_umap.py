from archimedes.functions.dataflow import output_path, input_path, input_var, json
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.plotting import px, pio, colors

scdata = sc.read(input_path('umap_anndata.h5ad'))

##### OUTPUT
fig = px.scatter(scdata.obsm['X_umap'], x=0, y=1, color_discrete_sequence=colors)

with open(output_path('umap.plotly.json'), 'w') as output_file:
    json.dump(json.loads(pio.to_json(fig)), output_file)
