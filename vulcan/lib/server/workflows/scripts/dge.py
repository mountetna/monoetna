from archimedes.functions.dataflow.DataIO import LocalData
from archimedes.functions.dataflow import input_path, input_var
from archimedes.functions.scanpy import scanpy as sc

def dge_on_features() -> Dict:
    rank_genes_groups_kwargs = {
        input_var('dge_method')
    }
    return rank_genes_groups_kwargs

def main(rank_genes_groups_kwargs: Dict) -> None:
    local_anndata = LocalData(input_path('leiden_anndata.h5ad'))
    local_anndata.io_wrapper(
        inner_func=sc.tl.rank_genes_groups,
        output_name='',
        **rank_genes_groups_kwargs
    )
