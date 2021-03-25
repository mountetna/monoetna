from archimedes.functions.dataflow import output_json, input_json
from archimedes.functions.magma import connect, question
from archimedes.functions.list import unique, flatten


seq_model_name = 'sc_seq'
seq_pool_model_name = 'sc_seq_pool'
magma = connect()

experiments = input_json('experiments')
tissues = input_json('tissues')

# Experiment and Tissue (AND logic)
filters = []
if len(experiments) > 0:
    filters.append( ['biospecimen_group', 'experiment', 'alias', '::in', experiments])
#### NEED TO TEST TISSUES BETTER ONCE ADDED
if len(tissues) > 0:
    filters.append( ['biospecimen_group', 'biospecimen_type', '::in', tissues] )

tube_records = unique(question(
    magma,
    [
        seq_model_name,
        [ '::has', 'raw_counts_h5'],
        *filters,
        '::all', '::identifier'
    ]))

output_json(tube_records, 'tube_recs')
