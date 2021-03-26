from archimedes.functions.dataflow import input_var, output_json
from archimedes.functions.magma import question, connect
from archimedes.functions.list import unique
from archimedes.functions.genes import GRCm38_ensembl93

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

fractions = question(
    magma,
    [
        seq_model_name,
        [ '::has', 'raw_counts_h5'],
        '::all',
        'cell_fraction'
    ]
)

color_options = {
    'Cluster': None,
    'Experiment': None,
    'Tissue': None,
    'Pool': None,
    'Biospecimen Group': None,
    'Gene': dict([ [ j[1], None ] for j in GRCm38_ensembl93 ])
}

output_json(color_options, 'color_options')
output_json(unique(experiments), 'experiments')
output_json(unique(tissues), 'tissues')
output_json(unique(fractions), 'fractions')
