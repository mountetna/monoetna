from archimedes.runner import RunRequest, StorageFile
import os.path

def get_fixture_path(fixture_path):
    return os.path.join(os.path.dirname(__file__), "fixtures", fixture_path)

def test_runner():
    request = RunRequest(
        script=open(get_fixture_path('test_runner_a.py'), 'r').read(),
        input_files=[StorageFile(host_path=get_fixture_path())],
    )