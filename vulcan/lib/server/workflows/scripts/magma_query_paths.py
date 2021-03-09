from archimedes.functions.dataflow import output_json, input_var
from archimedes.functions.magby import Magby
from archimedes.functions.environment import token

input_records = input_var('record_ids').split(", ")

# Improve these with magby.getModels or similar in future
seq_pool_model_name = 'sc_seq_pool'
seq_model_name = 'sc_seq'
h5_attr_name = 'raw_counts_h5'
project_name = "xcrs1" # <-- should be some sort of environment variable in future
magma = Magby.Magby(url='https://magma.ucsf.edu', token=token) # <-- <-- same as above comment

pool_records = [ rec for rec in input_records if rec.find('POOL') != -1 ]
tube_records = [ rec for rec in input_records if rec.find('POOL') == -1 ]

h5_locations = magma.query( projectName=project_name, queryTerms=[ seq_model_name,
    [ '::or',
       [ seq_pool_model_name, '::identifier', '::in', pool_records ],
       [ '::identifier', '::in', tube_records ]
    ], '::all', h5_attr_name, '::url'])
output_json(h5_locations, 'h5_locations')