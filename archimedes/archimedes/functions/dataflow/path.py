# Do not let these imports become usable within scripts themselves
import os as _os
import os.path as _os_path
from pandas import DataFrame, read_csv
import json


def input_path(name, inputs_env=_os.environ, inputs_dir=None):
    if inputs_dir is None:
        inputs_dir = inputs_env["INPUTS_DIR"]

    path = _os_path.join(inputs_dir, name)

    if not _os_path.exists:
        raise ValueError(f"Input key {name} does not exist for this cell")

    return path

def input_var(name, inputs_env=_os.environ, inputs_dir=None):
    return open(input_path(name, inputs_env, inputs_dir), 'r').read()

def input_tsv(name, inputs_env=_os.environ, inputs_dir=None):
    return read_csv( input_path(name, inputs_env, inputs_dir) , sep='\t')

def input_json(name, inputs_env=_os.environ, inputs_dir=None):
    return json.load(open(input_path(name, inputs_env, inputs_dir)))

def input_bool(name, inputs_env=_os.environ, inputs_dir=None):
    def str2bool(str):
        return str.lower() in ("yes", "true", "t", "1")
    return  str2bool(input_var(name, inputs_env, inputs_dir))

def output_path(name, outputs_env=_os.environ, outputs_dir=None):
    if outputs_dir is None:
        outputs_dir = outputs_env["OUTPUTS_DIR"]

    path = _os_path.join(outputs_dir, name)
    if outputs_env.get("ENFORCE_OUTPUTS_EXIST"):
        if not _os_path.exists(path):
            raise ValueError(f"Output key {name} does not exist for this cell.")

    return path

def output_tsv(data, name, outputs_env=_os.environ, outputs_dir=None):
    return data.to_csv( output_path(name, outputs_env, outputs_dir) , sep='\t')

def output_json(data, name, outputs_env=_os.environ, outputs_dir=None):
    with open(output_path(name, outputs_env, outputs_dir), 'w') as output_file:
        output_file.write(json.dumps(data))

def output_var(data, name, outputs_env=_os.environ, outputs_dir=None):
    with open(output_path(name, outputs_env, outputs_dir), 'w') as output_file:
        output_file.write(data)