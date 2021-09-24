from archimedes.functions.dataflow import output_path
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read('test_data/umap_workflow_anndata.h5ad')
scdata.write(output_path('scdata.h5ad'))