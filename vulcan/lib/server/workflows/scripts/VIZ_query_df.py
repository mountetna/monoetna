from archimedes.functions.dataflow import output_path, input_json
from archimedes.functions.utils import pandas as pd
from archimedes.functions.magma import connect, query_tsv
from archimedes.functions.environment import project_name

queryTerms = input_json('queryTerms')

magma = connect()
tsv_result = query_tsv(magma, project_name, queryTerms)

print(tsv_result)

tsv_result.to_json(output_path("data_frame"))
