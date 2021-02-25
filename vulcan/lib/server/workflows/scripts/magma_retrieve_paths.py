from archimedes.functions.dataflow import output_tsv, input_var
from archimedes.functions.magby import Magby

input_records = input_var('record_ids').split(", ")
token = input_var('token')

# Improve these with magby.getModels or similar in future
seq_pool_model_name = 'sc_rna_seq_pool' 
seq_model_name = 'sc_rna_seq'
h5_attr_name = ['raw_counts_h5']
project_name = "mvir1" # <-- should be some sort of environment variable in future
magma = Magby.Magby(url='https://magma.ucsf.edu', token=token) # <-- <-- same as above comment

singlular_record_ids = []
# Loop per record
for rec in input_records:
    
    # if a POOL, retrieve and add individual records
    if rec.find('POOL') != -1:
        ids_ret = magma.retrieve(
            projectName= project_name,
            modelName= seq_pool_model_name,
            recordNames= '[' + rec + ']',
            attributeNames= seq_model_name,
            dataType= 'meta')
        ids_str = ids_ret[seq_model_name]
        # split by ", "
        new_recs = ids_str.split(", ")
    else:
        # else, add id directly
        new_recs = rec

    singlular_record_ids.append(new_recs)

# Ensure uniqueness
singlular_record_ids = list(set(singlular_record_ids))

# Retrieve dataframe of record_id, raw_counts_h5 file location
h5_locations = magma.retrieve(
    projectName= project_name,
    modelName= seq_model_name,
    recordNames= singlular_record_ids,
    attributeNames= h5_attr_name,
    dataType= 'meta')

##### OUTPUT
output_tsv(h5_locations, 'h5_locations')