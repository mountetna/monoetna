from archimedes.functions.dataflow.DataIO import LocalData
from archimedes.functions.dataflow import input_path, input_var, input_json
from archimedes.functions.scanpy import scanpy as sc

def main() -> None:
    sc_rank_genes_groups_kwargs = input_json('selected_covariates')
    local_anndata = LocalData(input_path('leiden_anndata.h5ad'))
    local_anndata.io_wrapper(
        inner_func=sc.tl.rank_genes_groups,
        output_name='',
        **sc_rank_genes_groups_kwargs
    )
