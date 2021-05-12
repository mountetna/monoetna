from archimedes.functions.dataflow import output_path, input_path, output_json
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('leiden_anndata.h5ad'))

sc.tl.rank_genes_groups(scdata, 'leiden', method='wilcoxon')

top5 = pd.DataFrame(scdata.uns['rank_genes_groups']['names']).head(5)
output_json(top5.to_dict(), 'top5.json')