from archimedes.functions.dataflow import output_path, input_json, curl_data, tempfile, _os_path
from archimedes.functions.scanpy import scanpy as sc


def h5_dl_and_scanpy_import(tube_name, raw_counts_h5_mpath):
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmp_path = _os_path.join(tmpdirname,'new.h5')
        with open(tmp_path, 'wb') as tmp_file:
            tmp_file.write(curl_data(raw_counts_h5_mpath).content)
        adata = sc.read_10x_h5(tmp_path, gex_only = True)

    # add metadata named Library where it always equals tube_name
    adata.var_names_make_unique()
    adata.obs['Record_ID'] = tube_name

    return adata

data_tube_url = input_json('h5_locations')

# Initialize merged data, then loop through
merged_data = h5_dl_and_scanpy_import(data_tube_url[0][0], data_tube_url[0][1])

if len(data_tube_url)>1:
    # Loop per record
    for (tube_name, raw_counts_h5_mpath) in data_tube_url[1:]:
        # Convert from record name and metis path to the actual data
        new_data = h5_dl_and_scanpy_import(tube_name, raw_counts_h5_mpath)

        # Merge objects
        merged_data = merged_data.concatenate(new_data)

##### OUTPUT
merged_data.write(output_path('merged_anndata.h5ad'))
