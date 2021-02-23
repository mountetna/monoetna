from archimedes.functions.dataflow import output_path, input_path
from archimedes import scanpy

pools = input_var('pools')
# Convert these too
min_nCounts = open(input_path('min_nCounts'), 'r').read()
max_nCounts = open(input_path('max_nCounts'), 'r').read()
min_nFeatures = open(input_path('min_nFeatures'), 'r').read()
max_per_mito = open(input_path('max_per_mito'), 'r').read()
max_per_ribo = open(input_path('max_per_ribo'), 'r').read()

#### Probably outside of here, in R
## Read in .Rdata object and convert to h5ad OR export raw data as "10X h5" & metadata as tsv

##### magby retrieve gene_counts files....
## Obtain gene expression data + processed metadata if that exists.
## Fleshing this part out more fully will involve a discussion with Arjun about what's in the data library

# Assume output of magma querying will be aa set of raw_counts.h5 & metadata.tsv in some format for each record

merged_data = None
# Loop per record
for (tube_name, raw_counts_h5, metadata_file) in magma_results.items():
    # Convert from record name to the location of its h5 file

    ## scanpy processing
    # Read in
    adata = scanpy.read_10x_h5( # probably want a different scanpy read-in function
        raw_counts_h5,  # the directory with the `.mtx` file
        var_names='gene_ids'                # use gene symbols for the variable names (variables-axis index)
        )                              # write a cache file for faster subsequent reading

    ### Wrap in if (metadata not already available) statement
    metadata = DataFrame.from_tsv(metadata_file, sep='\t')
    # add metadata named Sample where that comes from the metadata
    adata.obs['sample'] = metadata['Sample.By.ABs'] # Figure out what it needs to be for this data, or expose as an input to the cell.
    
    # add metadata named Library where it always equals tube_name

    # Merge objects
    merged_data = merged_data.merge(adata) if merged_data else adata

##### OUTPUT
merged_data.write(output_path('merged_anndata'))

# Calculate mito QC metrics
merged_data.var['mt'] = merged_data.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
merged_data.var['rb'] = merged_data.var_names.str.startswith('RPS-|RPL-')  # annotate the group of ribosomal genes as 'rb' ## not sure if the syntax works.
scanpy.pp.calculate_qc_metrics(merged_data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
scanpy.pp.calculate_qc_metrics(merged_data, qc_vars=['rb'], percent_top=None, log1p=False, inplace=True)

# Perform all filtering
scanpy.pp.filter_cells(merged_data, min_genes=min_nFeatures, min_counts=min_nCounts, max_counts=max_nCounts)
merged_data = merged_data[merged_data.obs.pct_counts_mt < max_per_mito, :]
merged_data = merged_data[merged_data.obs.pct_counts_rb < max_per_ribo, :]

# Normalize: Counts depth then log transformation
scanpy.pp.normalize_total(merged_data, target_sum=1e4)
scanpy.pp.log1p(merged_data)

# Calculate genes to call highly variable
scanpy.pp.highly_variable_genes(merged_data, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Store a separate "assay" before we trim and scale 
merged_data.raw = merged_data

# Trim to variable genes
merged_data = merged_data[:, merged_data.var.highly_variable]

##### OUTPUT
merged_data.write(output_path('sc_rna_seq_normalized_data'))