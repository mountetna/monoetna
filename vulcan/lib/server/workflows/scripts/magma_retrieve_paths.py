from archimedes.functions.dataflow import output_tsv, input_var
from archimedes.functions.magby import Magby
from archimedes.functions.environment import token

input_records = input_var('record_ids').split(", ")

# Improve these with magby.getModels or similar in future
seq_pool_model_name = 'sc_rna_seq_pool' 
seq_model_name = 'sc_rna_seq'
h5_attr_name = 'raw_counts_h5'
project_name = "mvir1" # <-- should be some sort of environment variable in future
magma = Magby.Magby(url='https://magma.ucsf.edu', token=token) # <-- <-- same as above comment

singular_record_ids = []
# Loop per record
for rec in input_records:
    
    # if a POOL, retrieve and add individual records
    if rec.find('POOL') != -1:
        ids_ret = magma.retrieve(
            projectName= project_name,
            modelName= seq_pool_model_name,
            recordNames= rec,
            attributeNames= seq_model_name,
            dataType= 'meta')
        ids_str = ids_ret[seq_model_name]
        # split by ", "
        new_recs = ids_str[0].split(", ")

        singular_record_ids.extend(new_recs)
    else:
        # else, add id directly
        new_recs = rec

        singular_record_ids.append(new_recs)

# Ensure uniqueness
singular_record_ids = list(set(singular_record_ids))

# Retrieve dataframe of record_id, raw_counts_h5 file location
h5_locations = magma.retrieve(
    projectName= project_name,
    modelName= seq_model_name,
    recordNames= singular_record_ids,
    attributeNames= h5_attr_name,
    dataType= 'meta')

##### OUTPUT
output_tsv(h5_locations, 'h5_locations')