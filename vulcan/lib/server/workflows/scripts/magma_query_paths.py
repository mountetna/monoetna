from archimedes.functions.dataflow import output_json, input_var, input_json
from archimedes.functions.magma import connect, question
from archimedes.functions.environment import token, magma_host

input_records = input_json('record_ids')

# Improve these with magby.getModels or similar in future
seq_pool_model_name = 'sc_seq_pool'
seq_model_name = 'sc_seq'
h5_attr_name = 'raw_counts_h5'
project_name = "xcrs1"  # <-- should be some sort of environment variable in future
# <-- <-- same as above comment
magma = connect()

pool_records = [rec for rec in input_records if rec.find('POOL') != -1]
tube_records = [rec for rec in input_records if rec.find('POOL') == -1]

h5_locations = question(magma, [
    seq_model_name,
    [
        '::or',
        [ seq_pool_model_name, '::identifier', '::in', pool_records ],
        ['::identifier', '::in', tube_records]
    ],
    '::all', h5_attr_name, '::url'
], strip_identifiers=False)
output_json(h5_locations, 'h5_locations')
