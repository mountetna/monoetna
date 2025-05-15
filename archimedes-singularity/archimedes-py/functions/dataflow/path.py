# Do not let these imports become usable within scripts themselves
import os as _os
import os.path as _os_path
from pandas import DataFrame, read_csv
import json
from snakemake.script import snakemake

def input_path(key):
    if isinstance(key, str) and not key in snakemake.input:
        raise ValueError(f"The input key {key} does not exist for this job.")
    path = snakemake.input[key]
    if not _os_path.exists(path):
        raise ValueError(f"The input file for key {key} does not exist for this job.")
    return path

def input_var(key):
    return open(input_path(key), 'r').read()

def input_tsv(key):
    return read_csv(input_path(key) , sep='\t')

def input_json(key):
    return json.load(open(input_path(key)))

def input_bool(key):
    def str2bool(str):
        return str.lower() in ("yes", "y", "true", "t", "1")
    return  str2bool(input_var(key))

def param(key):
    if isinstance(key, str) and not key in snakemake.param:
        raise ValueError(f"The params key {key} does not exist for this job.")
    return snakemake@param[key]

def param_path(key):
    path = param(key)
    if not _os_path.exists(path):
        raise ValueError(f"No file at param key {key} path of {path}.")
    return path

def output_path(key):
    if isinstance(key, str) and not key in snakemake.output:
        raise ValueError(f"The output key {key} does not exist for this job.")
    return snakemake.output[key]

def output_tsv(data, key):
    return data.to_csv( output_path(key) , sep='\t')

def output_json(data, key):
    with open(output_path(key), 'w') as output_file:
        output_file.write(json.dumps(data))

def output_var(data, key):
    with open(output_path(key), 'w') as output_file:
        output_file.write(data)