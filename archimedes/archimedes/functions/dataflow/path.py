# Do not let these imports become usable within scripts themselves
import os as _os
import os.path as _os_path
from pandas import DataFrame


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
    return DataFrame.from_tsv( input_path(name, inputs_env, inputs_dir) , sep='\t')

def output_path(name, outputs_env=_os.environ, outputs_dir=None):
    if outputs_dir is None:
        outputs_dir = outputs_env["OUTPUTS_DIR"]

    path = _os_path.join(outputs_dir, name)
    if outputs_env.get("ENFORCE_OUTPUTS_EXIST"):
        if not _os_path.exists(path):
            raise ValueError(f"Output key {name} does not exist for this cell.")

    return path
