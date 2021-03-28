from archimedes.functions.dataflow import input_var, output_json
from archimedes.functions.magma import question, connect
from archimedes.functions.list import unique

seq_model_name = "sc_seq"
seq_pool_model_name = "sc_seq_pool"

magma = connect()

experiments = question(
    magma,
    [
        seq_model_name,
        [ '::has', 'raw_counts_h5'],
        '::all',
        'biospecimen_group', 'experiment', 'alias'
    ]
)

tissues = question(
    magma,
    [
        seq_model_name,
        [ '::has', 'raw_counts_h5'],
        '::all',
        'biospecimen_group', 'biospecimen_type'
    ]
)

output_json(unique(experiments), 'experiments')
output_json(unique(tissues), 'tissues')
