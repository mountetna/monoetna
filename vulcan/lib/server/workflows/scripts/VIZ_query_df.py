from archimedes.functions.dataflow import output_path, input_json, input_var, input_bool
from archimedes.functions.utils import json
from archimedes.functions.magma import connect, query_tsv
from archimedes.functions.environment import project_name

def wrap_brackets(str):
    if not (str.startswith("[") and str.endswith("]")):
        str = "[" + str + "]"
    return str

queryTerms = input_json('queryTerms')
user_columns = json.loads(wrap_brackets(input_var('user_columns')))
expand_matrices = input_bool('expand_matrices')

magma = connect()
tsv_result = query_tsv(magma, project_name,
    queryTerms=queryTerms,
    user_columns=user_columns,
    expand_matrices=expand_matrices,
    transpose=False)

tsv_result.to_json(output_path("data_frame"))
