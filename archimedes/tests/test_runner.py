from contextlib import contextmanager

from archimedes.runner import RunRequest, StorageFile, run, DockerIsolator, LocalIsolator
import os.path
import os

def read_resource(resource_path) -> str:
    return open(os.path.join(os.path.dirname(__file__), resource_path), 'r').read()

def host_storage_file(logical_name, app_path) -> StorageFile:
    return StorageFile(
        host_path=os.path.normpath(os.path.join(os.environ['HOST_DIR'], '..', 'archimedes', app_path)),
        logical_name=logical_name
    )

def local_storage_file(logical_name, app_path) -> StorageFile:
    return StorageFile(
        host_path=os.path.normpath(os.path.join(os.path.dirname(__file__), '..', app_path)),
        logical_name=logical_name
    )


@contextmanager
def expected_output(name):
    path = os.path.join(os.path.dirname(__file__), 'outputs', name)
    with open(path, 'w') as f:
        f.write('')

    yield f"tests/outputs/{name}", open(path, 'r')

    with open(path, 'w') as f:
        f.write('')

def test_runner():
    with expected_output('output_one') as (out_path, out_file):
        request = RunRequest(
            script=read_resource('fixtures/test_runner_a.py.donotrun'),
            input_files=[
                host_storage_file('input_one', 'tests/fixtures/runner_inputs/b')
            ],
            output_files=[
                host_storage_file('output_one', out_path)
            ],
        )

        result = run(request, DockerIsolator(), 4, False)
        assert result.error is None
        assert result.status == 'done'

        assert out_file.read() == "Some other input with more stuff in it"

    with expected_output('output_one') as (out_path, out_file):
        request = RunRequest(
            script=read_resource('fixtures/test_runner_a.py.donotrun'),
            input_files=[
                local_storage_file('input_one', 'tests/fixtures/runner_inputs/b')
            ],
            output_files=[
                local_storage_file('output_one', out_path)
            ],
        )

        result = run(request, LocalIsolator(), 4, False)
        assert result.error is None
        assert result.status == 'done'

        assert out_file.read() == "Some other input with more stuff in it"

    with expected_output('output_one') as (out_path, out_file):
        request = RunRequest(
            script=read_resource('fixtures/test_runner_a.py.donotrun'),
            input_files=[
                host_storage_file('input_one', 'tests/fixtures/runner_inputs/b')
            ],
            output_files=[
                host_storage_file('output_one', out_path)
            ]
        )

        result = run(request, DockerIsolator(), 4, False)
        assert result.error is None
        assert result.status == 'done'

        assert out_file.read() == "Some other input with more stuff in it"
