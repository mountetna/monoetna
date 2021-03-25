from archimedes.functions.dataflow import input_var, output_json
from archimedes.functions.environment import token, magma_host, project_name, app_env
from archimedes.functions.magby import Magby, query_extract
from archimedes.functions.list import unique, flatten

seq_model_name = "sc_seq"
seq_pool_model_name = "sc_seq_pool"
magma = Magby.Magby(url=magma_host, token=token,
                    verify=(app_env == 'production'))

experiments = magma.query(
    project_name,
    queryTerms=[
        "experiment",
        ["biospecimen_group", "::all",
        seq_model_name, "::all", '::has', 'raw_counts_h5'],
        "::all", 'alias'])
# flatten/unique
experiments = query_extract(experiments, 'alias')
experiments = ["No Selection"] + unique(flatten(experiments))

tissues = magma.query(
    project_name,
    queryTerms=[
        "biospecimen_group",
        [seq_model_name, "::all", '::has', 'raw_counts_h5'],
        "::all", 'biospecimen_type'])
tissues = query_extract(tissues, 'biospecimen_type')
tissues = ["No Selection"] + unique(flatten(tissues))

all_tubes = magma.query(
    project_name,
    queryTerms=[
        seq_model_name,
        ['::has', 'raw_counts_h5'],
        "::all", '::identifier']
        )['answer']
all_tubes = unique(flatten(all_tubes))


color_options = {
    'Experiment': {
        'High-fat': None,
        'Low-fat': None
    },
    'Tissue': {
        'Blood': None,
        'Tissue': None
    },
    'Gene': {
        'Gene1': None,
        'Gene2': None,
        'Gene3': None
    }
}

output_json(color_options, 'color_options')
output_json(experiments, 'experiments')
output_json(tissues, 'tissues')
output_json(all_tubes, 'all_tubes')
