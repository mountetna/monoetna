# Do not let these imports become usable within scripts themselves
import os as _os
import os.path as _os_path
from pandas import DataFrame, read_csv
import json

def input_path(key):

    # Todo: If key is string, check if the key exists, and error if not
    
    path = snakemake.input[key]

    if not _os_path.exists(path):
        raise ValueError(f"The input file for key {key} does not exist for this job.")
    
    return path

def input_var(key):
    return open(input_path(key), 'r').read()

def input_tsv(key):
    return read_csv( input_path(key) , sep='\t')

def input_json(key):
    return json.load(open(input_path(key)))

def input_bool(key):
    def str2bool(str):
        return str.lower() in ("yes", "y", "true", "t", "1")
    return  str2bool(input_var(key))

def output_path(key):

    # Todo: If key is string, check if the key exists, and error if not

    return snakemake.output[key]

def output_tsv(data, key):
    return data.to_csv( output_path(key) , sep='\t')

def output_json(data, key):
    with open(output_path(key), 'w') as output_file:
        output_file.write(json.dumps(data))

def output_var(data, key):
    with open(output_path(key), 'w') as output_file:
        output_file.write(data)