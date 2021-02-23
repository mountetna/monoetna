from archimedes.functions.dataflow import output_path, input_path
from archimedes.functions.metis import metis_file_path
from archimedes import scanpy

magma_output = input_tsv('magma_output')

merged_data = None
# Loop per record
for (tube_name, raw_counts_h5_mpath, metadata_file_mpath) in magma_output():
    # Convert from record name to the location of its h5 file

    ## scanpy processing
    # Read in
    adata = scanpy.read_10x_h5( # probably want a different scanpy read-in function
        metis_file_path(raw_counts_h5_mpath),  # the directory with the `.mtx` file
        var_names='gene_ids'                # use gene symbols for the variable names (variables-axis index)
        )                              # write a cache file for faster subsequent reading

    ### Wrap in if (metadata not already available) statement
    metadata = DataFrame.from_tsv(metis_file_path(metadata_file_mpath), sep='\t')
    # add metadata named Sample where that comes from the metadata
    adata.obs['sample'] = metadata['Sample.By.ABs'] # Figure out what it needs to be for this data, or expose as an input to the cell.
    
    # add metadata named Library where it always equals tube_name

    # Merge objects
    merged_data = merged_data.merge(adata) if merged_data else adata

##### OUTPUT
merged_data.write(output_path('merged_anndata'))