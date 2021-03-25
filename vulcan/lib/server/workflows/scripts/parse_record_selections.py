from archimedes.functions.dataflow import output_json, input_json
from archimedes.functions.magby import Magby
from archimedes.functions.environment import token, magma_host, app_env, project_name
from archimedes.functions.list import unique, flatten


seq_model_name = "sc_seq"
seq_pool_model_name = "sc_seq_pool"
magma = Magby.Magby(url=magma_host, token=token,
                    verify=(app_env == 'production'))

def defaultOnlyIfThere(x, default = "No Selection"):
    if default in set(x):
        x = [default]
    return x
experiments = defaultOnlyIfThere(input_json('experiments'))
tissues = defaultOnlyIfThere(input_json('tissues'))

# Experiment and Tissue (AND logic)
selection_terms = []
if experiments!=['No Selection']:
    selection_terms.append(
        ["biospecimen_group", "experiment",
        "alias", '::in', experiments])
#### NEED TO TEST TISSUES BETTER ONCE ADDED
if tissues!=['No Selection']:
    selection_terms.append(
        ["biospecimen_group", "biospecimen_type", '::in', tissues])

if len(selection_terms)!=0:
    if len(selection_terms)==1:
        selection_terms = selection_terms[0]
    else:
        selection_terms = ['::and'] + selection_terms
    
    tube_query = magma.query(
        projectName=project_name, 
        queryTerms=[
            seq_model_name,
            selection_terms,
            "::all", "::identifier"
        ])['answer']
    tube_records = unique(flatten(tube_query))
else:
    tube_records = input_json('all_tubes')

output_json(tube_records, 'tube_recs')
