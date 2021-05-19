from archimedes.functions.dataflow import input_path, input_json, output_json, parseModelAttr, buildTargetPath, curl_data, tempfile, _os_path
from archimedes.functions.magma import question, connect
from archimedes.functions.list import unique, flatten
from archimedes.functions.environment import project_name

# Retrieve scdata
# FILLER until actual data passing method is worked out
scdata = sc.read(input_path('umap_session_info'))
scdata.write(output_path('scdata'))
# with tempfile.TemporaryDirectory() as tmpdirname:
#     tmp_path = _os_path.join(tmpdirname,'scdata.h5ad')
#     with open(tmp_path, 'wb') as tmp_file:
#         tmp_file.write(curl_data(raw_counts_h5_mpath).content)
#     scdata = sc.read(tmp_path)

# Retrieve potential subset attributes
pdat = input_json("project_data")[project_name]
selection_atts = list(pdat['color_options'].keys()) + ['Cluster', 'Tube']
output_json(selection_atts, 'selection_atts')
