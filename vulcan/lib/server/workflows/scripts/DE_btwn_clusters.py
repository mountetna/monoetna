from archimedes.functions.dataflow import input_var, input_path, output_path, output_json
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.utils import pandas as pd
from archimedes.functions.list import unique

scdata = sc.read(input_path('leiden_anndata.h5ad'))

# Run differential expression between clusters
sc.tl.rank_genes_groups(scdata, 'leiden', method='wilcoxon')
DEdat = scdata.uns['rank_genes_groups']

# Output target "columns" of the dict output as csv
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

# Extract top 10 markers per cluster, ignoring given certain prefixes
names = DEdat['names']
out_len = 10

filters = input_var('ignore_prefixes')
if len(filters) > 0:
    def trim_prefix(list, prefix):
        for item in list[:]:
            if item.startswith(prefix):
                list.remove(item)
        return list
    def trim_prefixes(list, prefixes):
        for prefix in prefixes:
            list = trim_prefix(list, prefix)
        return list
    top10 = dict( [
        clust,
        trim_prefixes(list(names[clust]), filters.split(","))[:(out_len-1)]
    ] for clust in unique(scdata.obs['leiden']) )
else:
    top10 = pd.DataFrame(names).head(10).to_dict('list')

output_json(top10, 'top10.json')