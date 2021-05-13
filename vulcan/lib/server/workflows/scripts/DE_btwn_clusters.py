from archimedes.functions.dataflow import input_var, input_path, output_path, output_json
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.utils import pandas as pd
from archimedes.functions.list import unique

scdata = sc.read(input_path('leiden_anndata.h5ad'))
ignore_prefixes = input_var('ignore_prefixes')

sc.tl.rank_genes_groups(scdata, 'leiden', method='wilcoxon')

DEdat = scdata.uns['rank_genes_groups']

def DF_per_cluster(DEdat, cluster):
    DF = pd.DataFrame(
        dict([
            this_col,
            DEdat[this_col][cluster]
        ] for this_col in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
        )
    )
    DF['cluster'] = cluster
    return DF
pd.concat(
    (DF_per_cluster(DEdat, cluster)) for cluster in unique(scdata.obs['leiden'])
).to_csv(output_path('diffexp.csv'))

top10 = pd.DataFrame(DEdat['names']).head(10)
output_json(top10.to_dict(orient='list'), 'top10.json')