from archimedes.functions.dataflow import output_json, input_json
from archimedes.functions.magby import Magby
from archimedes.functions.environment import token, magma_host, app_env

# Expecting this may need to change...
experiments = input_json('experiments')
tissues = input_json('tissues')
pools = input_json('pools')
tubes = input_json('tubes')

### Need to make these actually do what they say!
def ignoreAllWhenOthers(x):
    return x
def removeBlank(x):
    return x
experiments = ignoreAllWhenOthers(experiments)
tissues = ignoreAllWhenOthers(tissues)
### May not need these actually?
pools = removeBlank(pools)
tubes = removeBlank(tubes)

# Setup query targets
expt_model = 'experiment'
expt_attr = '::' + 'identifier'
tiss_model = 'biospecimen'
tiss_attr = '::' + 'biospecimen_type'
project_name = "xcrs1"  # <-- should be some sort of environment variable in future

magma = Magby.Magby(url=magma_host, token=token,
                    verify=(app_env == 'production'))

# Build output
tube_records = []

# Experiment and Tissue (AND logic)
selection_terms = []
if experiments!='all':
    selection_terms.append([expt_model, expt_attr, '::in', experiments])
if tissues!='all':
    selection_terms.append([tiss_model, tiss_attr, '::in', tissues])

if len(selection_terms)!=0:
    if len(selection_terms)>1:
        selection_terms = ['::and', selection_terms]
    tube_records = magma.query(projectName=project_name, 
        queryTerms=[
            seq_model_name,
            selection_terms,
            '::all', '::identifier'])['answer']

# Pools and Tubes (OR logic / just add them!)
tube_records = tube_records.extend(
    magma.query(projectName=project_name, 
        queryTerms=[
            seq_model_name,
            ['::or',
                # Pools
                [seq_pool_model_name, '::identifier',
                '::in', pool_records],
                # Tubes
                ['::identifier',
                '::in', tube_records],
                '::all', '::identifier'])['answer']
    )

output_json(tube_records, 'tube_recs')
