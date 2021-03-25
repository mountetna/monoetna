from archimedes.functions.dataflow import input_var, output_json
from archimedes.functions.environment import token, magma_host, project_name, app_env
from archimedes.functions.magma import question, connect
from archimedes.functions.list import unique, flatten
from archimedes.functions.genes import GRCm38_ensembl93

def options(list, addDefault=True):
    return (['No Selection'] if addDefault else []) + unique(flatten(list))

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

all_tubes = question(
    magma,
    [
        seq_model_name,
        ['::has', 'raw_counts_h5'],
        "::all", '::identifier'
    ]
)


color_options = {
    'Cluster': None,
    'Experiment': None,
    'Tissue': None,
    'Pool': None,
    'Subject': None,
    'Gene': dict([ [ i, None ] for j in GRCm38_ensembl93 for i in j ])
}

output_json(color_options, 'color_options')
output_json(options(experiments), 'experiments')
output_json(options(tissues), 'tissues')
output_json(options(all_tubes, addDefault=False), 'all_tubes')
