from archimedes.functions.dataflow import output_path, input_path, output_json
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.utils import pandas as pd

scdata = sc.read(input_path('leiden_anndata.h5ad'))

sc.tl.rank_genes_groups(scdata, 'leiden', method='wilcoxon')

top10 = pd.DataFrame(scdata.uns['rank_genes_groups']['names']).head(10)
output_json(top10.to_dict(orient='list'), 'top10.json')