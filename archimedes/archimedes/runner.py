"""
This module implements a synchronous scheduler for running archimedes scripts in isolation.
In the future, the intention is for this to mostly disappear -- when we select a proper
scheduler, Vulcan will talk directly to that, and do so in an asynchronous fashion, submitting
jobs and reporting their status back to the frontend.

However, I believe this basic interface is stable enough to describe our use cases, so whatever
scheduler we move to should support a similar approach.  Perhaps this runner becomes just a library
that the scheduler itself uses locally to prepare the processes.  Or just a reference implementation
for moving thi logic elsewhere based on how we land on the scheduler question.
"""

import os.path
import os
import shutil
import signal
import subprocess
from dataclasses import dataclass
from dataclasses_json import dataclass_json
import time
from typing import List, ContextManager, TypeVar, Generic, Tuple, IO, Optional

from docker import DockerClient
from docker.types import Mount
from docker.models.containers import Container
from archimedes.checker import run as run_checker
import tempfile


@dataclass_json
@dataclass
class RunResult:
    status: str
    error: Optional[str]


@dataclass_json
@dataclass
class StorageFile:
    host_path: str
    logical_name: str


@dataclass_json
@dataclass
class RunRequest:
    script: str
    input_files: List[StorageFile]
    output_files: List[StorageFile]
    environment: List[str]
    image: str = "archimedes-base"
    isolator: str = "docker"


T = TypeVar('T')


class Isolator(Generic[T]):
    def start(self, request: RunRequest, exec_script_path) -> T:
        raise NotImplemented("Not implemented")

    def is_running(self, t: T) -> bool:
        raise NotImplemented("Not implemented")

    def stop(self, t: T, timeout=10):
        raise NotImplemented("Not implemented")

    def wait(self, t: T) -> int:
        raise NotImplemented("Not implemented")

    def get_stderr(self, t: T) -> str:
        raise NotImplemented("Not implemented")

    def reserve_exec_dir(self) -> ContextManager[str]:
        return tempfile.TemporaryDirectory()


class DockerIsolator(Isolator[Container]):
    docker_cli: DockerClient

    target_inputs_dir = "/inputs"
    target_outputs_dir = "/outputs"

    def __init__(self, docker_cli: DockerClient = DockerClient.from_env()):
        self.docker_cli = docker_cli

    def is_running(self, t: Container) -> bool:
        return t.status not in ['dead', 'exited']

    def stop(self, t: Container, timeout=10):
        return t.stop(timeout=timeout)

    def get_stderr(self, t: Container) -> str:
        return t.logs(stderr=True, stdout=False, tail='all')

    def wait(self, t: Container) -> int:
        result = t.wait()
        return result["StatusCode"]

    def start(self, request: RunRequest, exec_script_path: str) -> Container:
        mounts = [
                     Mount(target="/script.py", source=exec_script_path, read_only=True),
                 ] + [
                     Mount(target=os.path.join(self.target_inputs_dir, input_file.logical_name),
                           source=input_file.host_path)
                     for input_file in request.input_files
                 ] + [
                     Mount(target=os.path.join(self.target_outputs_dir, output_file.logical_name),
                           source=output_file.host_path)
                     for output_file in request.output_files
                 ]

        environment = request.environment + [
            f"INPUTS_DIR={self.target_inputs_dir}",
            f"OUTPUTS_DIR={self.target_outputs_dir}",
            f"ENFORCE_OUTPUTS_EXIST=1",
        ]

        # Could set options like cpu_quote, mem_limit, restrict network access,
        # etc, to further lockdown the task.
        return self.docker_cli.containers.run(
            request.image,
            f"poetry run /script.py",
            auto_remove=True,
            detach=True,
            environment=environment,
            # Arbitrary limit of 2G for any other buffered data.
            storage_opt=dict(size="2G"),
            mounts=mounts,
        )


LocalProcess = Tuple[subprocess.Popen, IO[bin], str]


class LocalIsolator(Isolator[LocalProcess]):
    root_dir: str

    def is_running(self, t: LocalProcess) -> bool:
        p, io, f = t
        return p.returncode is None

    def stop(self, t: LocalProcess, timeout=10):
        p, io, f = t

        done = time.time() + timeout
        while self.is_running(t) and time.time() < done:
            p.send_signal(signal.SIGTERM)
            time.sleep(1)

        p.kill()

    def get_stderr(self, t: LocalProcess) -> str:
        p, io, f = t
        return open(f, "rb").read().decode('utf-8')

    def wait(self, t: LocalProcess) -> int:
        p, io, f = t
        return p.wait()

    def start(self, request: RunRequest, exec_script_path: str) -> LocalProcess:
        inputs_dir = os.path.join(exec_script_path, "inputs")
        outputs_dir = os.path.join(exec_script_path, "outputs")
        stderr_log = os.path.join(exec_script_path, 'stderr')

        os.mkdir(inputs_dir)
        os.mkdir(outputs_dir)

        for input_file in request.input_files:
            os.symlink(input_file.host_path, os.path.join(inputs_dir, input_file.logical_name))

        for output_file in request.output_files:
            os.symlink(output_file.host_path, os.path.join(outputs_dir, output_file.logical_name))

        environment = request.environment + [
            f"INPUTS_DIR={inputs_dir}",
            f"OUTPUTS_DIR={outputs_dir}",
            f"ENFORCE_OUTPUTS_EXIST=1",
        ]

        stderr_log_io = open(stderr_log, 'wb')
        return subprocess.Popen(
            [shutil.which('python'), exec_script_path],
            stderr=stderr_log_io,
            env={k: v for k, v in [s.split('=', maxsplit=1) for s in environment]},
        ), stderr_log_io, stderr_log


def run(request: RunRequest, isolator: Isolator[T]) -> RunResult:
    res = RunResult(status='running', error=None)

    with isolator.reserve_exec_dir() as exec_dir:
        script_path = os.path.join(exec_dir, "script.py")
        with open(script_path, "w") as script_file:
            script_file.write(request.script)

        if not run_checker([script_path]):
            raise ValueError(
                "Script did not conform to the archimedes dsl requirements."
            )

        process = isolator.start(request, script_path)

        # Give cells up to 5 minutes.  In the future, with a proper async executor we would want to allow up to
        # a much larger amount of time.
        done = time.time() + 60 * 5

        while isolator.is_running(process) and time.time() < done:
            time.sleep(1)

        if isolator.is_running(process):
            isolator.stop(process)

        code = isolator.wait(process)

        if code != 0:
            res.status = 'failed'
            res.error = isolator.get_stderr(process)
        else:
            res.status = 'done'

    return res


def main():
    import sys

    request: RunRequest = RunRequest.schema().loads(sys.stdin.read())
    if request.isolator == 'docker':
        result = run(request, DockerIsolator())
        print(RunResult.schema().dumps(result))
