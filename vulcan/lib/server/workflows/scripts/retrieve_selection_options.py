from archimedes.functions.dataflow import input_var, output_json
from archimedes.functions.environment import token, magma_host, project_name, app_env
from archimedes.functions.magby import Magby
from archimedes.functions.list import unique, flatten

seq_model_name = "sc_seq"
seq_pool_model_name = "sc_seq_pool"
magma = Magby.Magby(url=magma_host, token=token,
                    verify=(app_env == 'production'))

experiments = magma.query(
    project_name,
    queryTerms=[
        "experiment",
        ["subject", "::all",
        "biospecimen", "::all",
        seq_model_name, "::all", '::has', 'raw_counts_h5'],
        "::all", '::identifier'])['answer']
# add "all" and flatten/unique
experiments = ["all"] + unique(flatten(experiments))

tissues = magma.query(
    project_name,
    queryTerms=[
        "biospecimen",
        [seq_model_name, "::all", '::has', 'raw_counts_h5'],
        "::all", 'biospecimen_type'])['answer']
tissues = ["all"] + unique(flatten(tissues))

pools = magma.query(
    project_name,
    queryTerms=[
        seq_pool_model_name,
        [seq_model_name, "::all",
        '::has', 'raw_counts_h5'],
        "::all", '::identifier'])['answer']
# Add "" and flatten/unique
pools = ["NO DIRECT SELECTIONS"] + unique(flatten(pools))

recs = magma.query(
    project_name,
    queryTerms=[
        seq_model_name,
        ['::has', 'raw_counts_h5'],
        "::all", '::identifier'])['answer']
recs = ["NO DIRECT SELECTIONS"] + unique(flatten(recs))

output_json(experiments, 'experiments')
output_json(tissues, 'tissues')
output_json(pools, 'pools')
output_json(recs, 'records')
