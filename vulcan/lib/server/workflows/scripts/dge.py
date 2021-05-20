from archimedes.functions.dataflow.DataIO import LocalData
from archimedes.functions.dataflow import input_path
from archimedes.functions.scanpy import scanpy as sc

def dge_on_features(scdata: AnnData, **kwargs) -> AnnData:
    pass

def main(rank_genes_groups_kwargs: Dict):
    local_anndata = LocalData(input_path('leiden_anndata.h5ad'))
    local_anndata.io_wrapper(
        inner_func=sc.tl.rank_genes_groups,
        output_name='',
        **rank_genes_groups_kwargs
    )
