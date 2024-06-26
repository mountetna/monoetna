from archimedes.functions.dataflow import input_var, input_path, output_path, output_json
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.utils import pandas as pd
from archimedes.functions.list import unique

scdata = sc.read(input_path('leiden_anndata.h5ad'))

# Run differential expression between clusters
dge_method = input_var('dge_method')
sc.tl.rank_genes_groups(scdata, groupby='leiden', method=dge_method, pts=True, corr_method='benjamini-hochberg')
DEdat = scdata.uns['rank_genes_groups']

# Output the anndata object
scdata.write(output_path('umap_workflow_anndata.h5ad'))

# Output target "columns" of the dict output as csv
def collect_DF_per_cluster(DEdat, cluster):
    keys = ['names', 'pvals', 'pvals_adj', 'logfoldchanges']
    DF = pd.DataFrame(
        dict(
            [
                this_col,
                DEdat[this_col][cluster]
            ] for this_col in keys
        )
    )
    ## Add pts (percent expressed), ensuring proper ordering
    names = list(DF['names'])
    for key in ['pts', 'pts_rest']:
        DF[key] = list(DEdat[key].reindex(names)[cluster])
    DF['cluster'] = cluster
    return DF
pd.concat(
    (collect_DF_per_cluster(DEdat, cluster)) for cluster in unique(scdata.obs['leiden'])
).to_csv(output_path('cluster_diffexp.csv'))

# Extract top 10 markers per cluster, ignoring given certain prefixes
names = DEdat['names']
out_len = 10

filters = input_var('ignore_prefixes').lower()
if len(filters) > 0:
    def trim_prefix(list, prefix):
        for item in list[:]:
            if item.lower().startswith(prefix):
                list.remove(item)
        return list
    def trim_prefixes(list, prefixes):
        for prefix in prefixes:
            list = trim_prefix(list, prefix)
        return list
    top10 = dict( [
        clust,
        trim_prefixes(list(names[clust]), filters.split(","))[:(out_len)]
    ] for clust in unique(scdata.obs['leiden']) )
else:
    top10 = pd.DataFrame(names).head(out_len).to_dict('list')

output_json(top10, 'top10.json')