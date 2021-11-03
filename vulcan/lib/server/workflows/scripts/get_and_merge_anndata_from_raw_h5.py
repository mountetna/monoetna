from archimedes.functions.dataflow import output_path, input_json, curl_data, tempfile, _os_path, buildTargetPath, parseModelAttr
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.magma import connect, question
from archimedes.functions.environment import token, magma_host, project_name
from archimedes.functions.utils import re
from archimedes.functions.list import unique

input_records = input_json('record_ids')

### Query paths to raw data for requested records
pdat = input_json("project_data")[project_name]
seq_target = parseModelAttr(pdat['seq_h5_counts_data'])

magma = connect()

data_tube_url = question(
    magma, 
    [
        seq_target['model'],
        ['::identifier', '::in', input_records],
        '::all', seq_target['attribute'], '::url'
    ],
    strip_identifiers=False)

### Initialize merged data, then loop through
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

merged_data = h5_dl_and_scanpy_import(data_tube_url[0][0], data_tube_url[0][1])

if len(data_tube_url)>1:
    # Loop per record
    for (tube_name, raw_counts_h5_mpath) in data_tube_url[1:]:
        # Convert from record name and metis path to the actual data
        new_data = h5_dl_and_scanpy_import(tube_name, raw_counts_h5_mpath)

        # Merge objects
        merged_data = merged_data.concatenate(new_data)

### Add metadata for all color-data
pdat = input_json("project_data")[project_name]
color_options = pdat['color_options']

# Prep magma querying
seq_model = parseModelAttr(pdat['seq_h5_counts_data'])['model']
magma = connect()
def get(ids, value):
    ids = ids if type(ids) == list else ids.tolist()
    values = dict(question(magma, [
        seq_model,
        [ '::identifier', '::in', ids ],
        '::all',
        *value
    ], strip_identifiers=False))

    return [ values.get(id, None) for id in ids ]

# Add data
recs = merged_data.obs[ 'Record_ID' ].values
for color_by in color_options.keys():
    rec_data = dict(map(
        lambda x: [x, get([x], buildTargetPath(color_options[color_by], pdat))[0]],
        unique(recs) ))
    label = re.sub(" ", "_", color_by)
    #raise Exception(rec_data)
    merged_data.obs[ label ] = list( map( lambda x: rec_data[x], recs) )

### OUTPUT
merged_data.write(output_path('merged_anndata.h5ad'))
