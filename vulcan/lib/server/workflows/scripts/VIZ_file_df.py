from archimedes.functions.dataflow import output_path, input_json
from archimedes.functions.utils import pandas
from archimedes.functions.environment import project_name, metis_host, token
from archimedes.functions.etna import Metis, File, TokenAuth
from archimedes.functions.utils import re

metis = Metis(TokenAuth(token), metis_host)
fileInfo = input_json('data_table_file')
file = File(
    file_path=fileInfo['path'], project_name=project_name, bucket_name=fileInfo['bucket'], download_url = None
)

separator = "\t"
if re.search("csv$", file.file_path) is not None:
    separator = ","

with metis.open_file(file) as open_file:
    df = pandas.read_table(open_file, sep=separator, engine="python")
df.to_json(output_path("data_frame"))
