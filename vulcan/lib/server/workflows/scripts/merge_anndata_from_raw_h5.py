from archimedes.functions.dataflow import output_path, input_tsv, curl_data
from archimedes.functions.scanpy import scanpy as sc

data_urls = input_tsv('h5_locations')

def scanpy_import(tube_name, raw_counts_h5_mpath):
    ## scanpy Read-in
    # Read in
    adata = sc.read_10x_h5(
        curl_data(raw_counts_h5_mpath)
        )

    # add metadata named Library where it always equals tube_name
    adata.obs['Record_ID'] = tube_name

    return adata

# Initialize merged data, then loop through
merged_data = scanpy_import(data_urls.tube_name[0], data_urls.raw_counts_h5[0])

if data_urls['tube_name'].length()>1:
    # Loop per record
    for (tube_name, raw_counts_h5_mpath) in data_urls[1:].iterrows():
        # Convert from record name and metis path to the actual data
        new_data = scanpy_import(tube_name, raw_counts_h5)

        # Merge objects
        merged_data = merged_data.merge(new_data) if merged_data else adata

##### OUTPUT
merged_data.write(output_path('merged_anndata'))