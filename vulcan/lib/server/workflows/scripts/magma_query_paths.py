from archimedes.functions.dataflow import output_json, input_var, input_json, parseModelAttr, buildTargetPath
from archimedes.functions.magma import connect, question
from archimedes.functions.environment import token, magma_host, project_name

input_records = input_json('record_ids')

pdat = input_json("project_data")[project_name]
seq_target = parseModelAttr(pdat['seq_h5_counts_data'])

magma = connect()

h5_locations = question(
    magma, 
    [
        seq_target['model'],
        ['::identifier', '::in', input_records],
        '::all', seq_target['attribute'], '::url'
    ],
    strip_identifiers=False)
output_json(h5_locations, 'h5_locations')
