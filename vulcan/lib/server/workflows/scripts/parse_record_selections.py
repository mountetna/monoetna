from archimedes.functions.dataflow import output_json, input_json
from archimedes.functions.magby import Magby
from archimedes.functions.environment import token, magma_host, app_env
from archimedes.functions.list import unique, flatten

# Expecting this may need to change...
experiments = input_json('experiments')
tissues = input_json('tissues')
pools = input_json('pools')
tubes = input_json('tubes')

seq_model_name = "sc_seq"
seq_pool_model_name = "sc_seq_pool"

def defaultOnlyIfThere(x, default = "No Selection"):
    if default in set(x):
        x = [default]
    return x
experiments = defaultOnlyIfThere(experiments)
tissues = defaultOnlyIfThere(tissues)
pools = defaultOnlyIfThere(pools)
tubes = defaultOnlyIfThere(tubes)

print("expts:", experiments)
print("tissues:", tissues)
print("pools:", pools)
print("tubes:", tubes)

# Setup query targets
# expt_model = 'experiment'
# expt_attr = '::' + 'identifier'
# tiss_model = 'biospecimen'
# tiss_attr = 'biospecimen_type'
project_name = "xcrs1"

magma = Magby.Magby(url=magma_host, token=token,
                    verify=(app_env == 'production'))

# Build output
tube_records = []

# Experiment and Tissue (AND logic)
selection_terms = []
if experiments!=['No Selection']:
    selection_terms.append(
        ["biospecimen", "subject", "experiment",
        "::identifier", '::in', experiments])
#### NEED TO TEST TISSUES ONCE ADDED
if tissues!=['No Selection']:
    selection_terms.append(
        ["biospecimen", "biospecimen_type", '::in', tissues])

if len(selection_terms)!=0:
    if len(selection_terms)>=2:
        selection_terms = ['::and'] + selection_terms + ['::all', '::identifier']
        # tube_records.extend(
        print("query 2+:",
            magma.query(
                projectName=project_name, 
                queryTerms=[
                    seq_model_name,
                    selection_terms
                ])#['answer']
        )
    else:
        tube_records.extend(
            magma.query(
                projectName=project_name, 
                queryTerms=[
                    seq_model_name,
                    selection_terms[0], # double-list to single
                    "::all", "::identifier"
                ])['answer']
        )

# Pools and Tubes (OR logic / just add them!)
selection_terms = []
#### NEED TO TEST POOLS ONCE ADDED
if pools != ['No Selection']:
    tube_records.extend(
        magma.query(projectName=project_name, 
            queryTerms=[
                seq_model_name,
                [seq_pool_model_name, '::identifier',
                        '::in', pools],
                '::all', '::identifier'
            ])['answer']
        )
if tubes != ['No Selection']:
    tube_records.extend(tubes)        

print("before trim:", tube_records)

tube_records = unique(flatten(tube_records))

print("out:", tube_records)

output_json(tube_records, 'tube_recs')
