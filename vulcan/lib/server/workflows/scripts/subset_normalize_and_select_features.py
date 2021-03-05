from archimedes.functions.dataflow import output_path, input_path, input_var
from archimedes.functions.scanpy import scanpy as sc

scdata = sc.read(input_path('merged_anndata.h5ad'))
# Convert these too
min_nCounts = int(input_var('min_nCounts'))
max_nCounts = int(input_var('max_nCounts'))
min_nFeatures = int(input_var('min_nFeatures'))
max_per_mito = int(input_var('max_per_mito'))
max_per_ribo = int(input_var('max_per_ribo'))

# Calculate mito QC metrics
scdata.var['mt'] = scdata.var_names.str.startswith('Mt-')  # annotate the group of mitochondrial genes as 'mt'
scdata.var['rb'] = scdata.var_names.str.startswith('Rps-|Rpl-')  # annotate the group of ribosomal genes as 'rb' ## not sure if the syntax works.
sc.pp.calculate_qc_metrics(scdata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.calculate_qc_metrics(scdata, qc_vars=['rb'], percent_top=None, log1p=False, inplace=True)

# Perform all filtering
sc.pp.filter_cells(scdata, min_genes=min_nFeatures)
sc.pp.filter_cells(scdata, min_counts=min_nCounts)
sc.pp.filter_cells(scdata, max_counts=max_nCounts)
scdata = scdata[scdata.obs.pct_counts_mt < max_per_mito, :]
scdata = scdata[scdata.obs.pct_counts_rb < max_per_ribo, :]

# Normalize: Counts depth then log transformation
sc.pp.normalize_total(scdata, target_sum=1e4)
sc.pp.log1p(scdata)

# Calculate genes to call highly variable
sc.pp.highly_variable_genes(scdata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Store a separate "assay" before we trim and scale 
scdata.raw = scdata

# Trim to variable genes
scdata = scdata[:, scdata.var.highly_variable]

##### OUTPUT
scdata.write(output_path('sc_rna_seq_normalized_data'))