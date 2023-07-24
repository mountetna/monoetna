from archimedes.functions.dataflow import output_path, input_json, input_var, input_bool
from archimedes.functions.utils import json, pandas
from archimedes.functions.magma import connect, query_tsv
from archimedes.functions.environment import project_name, metis_host, token
from archimedes.functions.etna import Metis, File, TokenAuth

def wrap_brackets(str):
    if not (str.startswith("[") and str.endswith("]")):
        str = "[" + str + "]"
    return str

useQuery = input_bool('useQuery')

if useQuery:
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
else:
    metis = Metis(TokenAuth(token), metis_host)
    fileInfo = input_json('fileInfo')
    file = File(
        file_path=fileInfo['path'], project_name=project_name, bucket_name=fileInfo['bucket'], download_url = None
    )
    file_reader = pandas.read_csv if fileInfo['path'].endswith('csv') else pandas.read_table
    with metis.open_file(file) as open_file:
        df = file_reader(open_file)
    df.to_json(output_path("data_frame"))
