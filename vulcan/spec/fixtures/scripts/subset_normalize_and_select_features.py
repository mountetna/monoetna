from archimedes.functions.dataflow import output_path, input_path
from archimedes import scanpy

pools = open(input_path('sc_rna_seq_pools'), 'r').read()
min_nCounts = open(input_path('min_nCounts'), 'r').read()
max_nCounts = open(input_path('max_nCounts'), 'r').read()
min_nFeatures = open(input_path('min_nFeatures'), 'r').read()
max_per_mito = open(input_path('max_per_mito'), 'r').read()
max_per_ribo = open(input_path('max_per_ribo'), 'r').read()

##### magby retrieve gene_counts files....
## Empty Stub


##### scanpy processing

# Read in
adata = scanpy.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_ids',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading

# Filter cells by min_genes
scanpy.pp.filter_cells(adata, min_genes=200)

# Calculate mito QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
scanpy.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter on max genes and mitochondria
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Normalize: Counts depth then log transformation
scanpy.pp.normalize_total(adata, target_sum=1e4)
scanpy.pp.log1p(adata)

# Calculate genes to call highly variable
scanpy.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Store a separate "assay" before we trim and scale 
adata.raw = adata

# Trim to variable genes
adata = adata[:, adata.var.highly_variable]

##### OUTPUT
adata.write(output_path('sc_rna_seq_normalized_data'))