# Do not let these imports become usable within scripts themselves
import os as _os
import os.path as _os_path


def open_input(name, inputs_env=_os.environ, inputs_dir=None):
    if inputs_dir is None:
        inputs_dir = inputs_env.get("INPUTS_DIR") or _os.getcwd()
        if inputs_dir is None:
            inputs_dir = _os_path.join(_os.getcwd(), "inputs")
            _os.mkdir(inputs_dir)

    return open(_os_path.join(inputs_dir, name), "rb")


def open_output(name, outputs_env=_os.environ, outputs_dir=None):
    if outputs_dir is None:
        outputs_dir = outputs_env.get("OUTPUTS_DIR")
        if outputs_dir is None:
            outputs_dir = _os_path.join(_os.getcwd(), "outputs")
            _os.mkdir(outputs_dir)

    path = _os_path.join(outputs_dir, name)
    if outputs_env.get("ENFORCE_OUTPUTS_EXIST"):
        if not _os_path.exists(path):
            raise ValueError(f"Output key {name} was not specified for this script.")

    return open(path, "wb")
