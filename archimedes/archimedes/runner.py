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
from dataclasses import dataclass, field
from dataclasses_json import dataclass_json
import time
from typing import List, ContextManager, TypeVar, Generic, Tuple, IO, Optional

from docker import DockerClient
from docker.errors import NotFound
from docker.types import Mount
from docker.models.containers import Container
from kubernetes import client, config, watch
from kubernetes.client.api import core_v1_api
from kubernetes.client.rest import ApiException
from archimedes.checker import run as run_checker
import tempfile
import sys
import secrets


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
    input_files: List[StorageFile] = field(default_factory=list)
    output_files: List[StorageFile] = field(default_factory=list)
    environment: List[str] = field(default_factory=list)
    image: str = "archimedes"
    isolator: str = "docker"
    interpreter: str = "python"

    def run_with(self, interp: str):
        return self.interpreter == interp

    @property
    def target_file(self):
        if self.interpreter == 'node':
            return "script.js"
        if self.interpreter == 'r':
            return "script.r"
        else:
            return "script.py"

    @property
    def target_path(self):
        return f"/{self.target_file}"


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

    # This needs to be a unique, top level directory to enable sane
    # volume directory sharing.
    exec_dir = '/archimedes-exec'
    target_inputs_dir = "/inputs"
    target_outputs_dir = "/outputs"

    def __init__(self, docker_cli: Optional[DockerClient] = None):
        self.docker_cli = docker_cli or DockerClient.from_env()

    def is_running(self, t: Container) -> bool:
        try:
            t.reload()
            return t.status == 'running'
        except NotFound:
            return False

    def stop(self, t: Container, timeout=10):
        try:
            t.reload()
            return t.stop(timeout=timeout)
        except NotFound:
            return

    def get_stderr(self, t: Container) -> str:
        try:
            t.reload()
            return t.logs(stderr=True, stdout=True, tail='all')
        except NotFound:
            return 'Process disappeared'

    def wait(self, t: Container) -> int:
        result = t.wait()
        return result["StatusCode"]

    def _is_self_container(self):
        proc_file = '/proc/self/cgroup'

        if os.path.exists('/.dockerenv'):
            return True

        if not os.path.exists(proc_file):
            return False

        return b"/docker" in open(proc_file, 'rb').read()

    def start(self, request: RunRequest, exec_script_path: str) -> Container:
        mounts =  [
                     Mount(target=os.path.join(self.target_inputs_dir, input_file.logical_name),
                           type='bind',
                           source=input_file.host_path)
                     for input_file in request.input_files
                 ] + [
                     Mount(target=os.path.join(self.target_outputs_dir, output_file.logical_name),
                           type='bind',
                           source=output_file.host_path)
                     for output_file in request.output_files
                 ]

        environment = request.environment + [
            f"INPUTS_DIR={self.target_inputs_dir}",
            f"OUTPUTS_DIR={self.target_outputs_dir}",
            f"ENFORCE_OUTPUTS_EXIST=1",
        ]

        volumes_from = []

        if not self._is_self_container():
            mounts += [
                Mount(
                    target=request.target_path,
                    type='bind',
                    source=exec_script_path,
                    read_only=True
                )
            ]
            script_path = request.target_path
        else:
            volumes_from += [ os.environ['HOSTNAME'] ]
            script_path = exec_script_path

        # Could set options like cpu_quote, mem_limit, restrict network access,
        # etc, to further lockdown the task.
        params = {
            "detach": True,
            "environment": environment,
            "volumes_from": volumes_from,
            "mounts": mounts,
            "network_mode": 'host',
        }

        return self.docker_cli.containers.run(
            request.image, self._cmd(request, script_path),
            **params
        )

    def reserve_exec_dir(self) -> ContextManager[str]:
        if not os.path.exists(self.exec_dir):
            os.mkdir(self.exec_dir)

        return tempfile.TemporaryDirectory(dir=self.exec_dir)

    def _cmd(self, request: RunRequest, exec_script_path: str):
        if request.run_with("node"):
            exec_dir = exec_script_path.replace(request.target_file, "")
            return ["/bin/bash", "-c", f"cd {exec_dir}; cp /app/package.json package.json; npm install; node {exec_script_path}"]
        elif request.run_with("r"):
            return f"Rscript {exec_script_path}"
        else:
            return f"poetry run python {exec_script_path}"

class KubernetesIsolator(Isolator[str]):
    # This needs to be a unique, top level directory to enable sane
    # volume directory sharing.
    exec_dir = '/archimedes-exec'
    target_inputs_dir = "/inputs"
    target_outputs_dir = "/outputs"
    result = 1

    def get_container_state(self, pod):
        container = pod.status.container_statuses[-1]
        return [k for k in container.state.to_dict().items() if k[1] is not None][0]

    def pod_container_status(self, t: str, namespace='default'):
        pod = None
        try:
            pod = self.kub_cli.read_namespaced_pod(t, namespace=namespace)
            return self.get_container_state(pod)
        except:
            return ("Unknown", pod)
    
    def delete_pod(self, t: str, namespace='default', timeout=60, maybe_repeat=False):
        try:
            self.kub_cli.delete_namespaced_pod(t, namespace=namespace, grace_period_seconds=timeout)
            return
        except:
            if maybe_repeat:
                print("Error deleting pod {t}, but might have been deleted already.", file=sys.stderr)
            raise Exception(f"Error deleting pod {t}.")

    def pod_and_log(self, t: str, namespace='default', pod_only: bool = False):
        try:
            log = ''
            if not pod_only:
                log = self.kub_cli.read_namespaced_pod_log(t, namespace=namespace)
            pod = self.kub_cli.read_namespaced_pod(t, namespace=namespace)
            return (pod, log)
        except:
            raise Exception(f"Error reading pod {t}.")

    def __init__(self, kub_cli: Optional[DockerClient] = None):
        if kub_cli is None:
            # ToDo: config specified location!
            config.load_kube_config('test_data/config')
            self.kub_cli = core_v1_api.CoreV1Api()
        else:
            self.kub_cli = kub_cli

    def is_running(self, t: str, namespace='default') -> bool:
        if self.pod_and_log(t, namespace=namespace, pod_only=True)[0].status.phase == "Pending":
            return True
        container_state = self.pod_container_status(t, namespace=namespace)
        if container_state[0] == "Unknown":
            raise Exception("Unknown container state", container_state[1])
        return container_state[0] in ["waiting", "running"]

    ## Because there is no analagous concept for the wait() goal of 'container is wrapping up, await that and report outcome'...
    # And because the usage is 'stop() on timeout, but always call wait() at end and use that to get final status code'...
    # stop(): will instead do nothing
    # wait(): will record the current / last container status, then perform the stop() equivalent
    def wait(self, t: str, timeout=60) -> Tuple:
        if self.pod_container_status(t)[0] == "terminated":
            pod,log=self.pod_and_log(t)
            return (pod.status.container_statuses[-1].state.terminated.exit_code, log)
    
        log = self.pod_and_log(t)[1]
        pod = self.delete_pod(t, grace_period_seconds=timeout)
        # await grace_period, hoping for completion
        done = time.time() + timeout
        while self.is_running(t) and time.time() < done:
            pod,log = self.pod_and_log(t)
            time.sleep(1)
        
        if self.get_container_state(t)[0] == "terminated":
            return (pod.status.container_statuses[-1].state.terminated.exit_code, log)
        
        return ("Pod killed", log)

    # Handled within wait()
    def stop(self, t: str, timeout=60):
        return

    # Handled within wait(), should never be called
    def get_stderr(self, t: str) -> str:
        raise NotImplemented("Not implemented")

    # def _is_self_container(self):
    #     proc_file = '/proc/self/cgroup'

    #     if os.path.exists('/.dockerenv'):
    #         return True

    #     if not os.path.exists(proc_file):
    #         return False

    #     return b"/docker" in open(proc_file, 'rb').read()

    def kubernetes_run(self, image, _cmd, container_params, pod_params):

        pod_name = 'archimedes-job-' + secrets.token_hex(10)

        pod_manifest = {
            'apiVersion': 'v1',
            'kind': 'Pod',
            'metadata': {
                'name': pod_name
            },
            'spec': {
                'restart_policy': 'Never',
                'image_pull_secrets':[{
                    'name': 'regcred'
                }],
                "hostNetwork": True,
                'containers': [{
                    'image': image,
                    'image_pull_policy': 'IfNotPresent',
                    'name': 'archimedes-job',
                    "args": _cmd,
                    **container_params
                }],
                **pod_params
            }
        }
        # Returns V1Pod
        resp =  self.kub_cli.create_namespaced_pod(body=pod_manifest,
                                                  namespace='default')
        time.sleep(1)
        return pod_name

    def start(self, request: RunRequest, exec_script_path: str):
        mounts = [
                    {
                        "name": f"input-{input_file.logical_name}",
                        "mount_path": os.path.join(self.target_inputs_dir, input_file.logical_name),
                        "host_path": input_file.host_path,
                        "read_only": True,
                    }
                    for input_file in request.input_files
                 ] + [
                        {
                            "name": f"output-{output_file.logical_name}",
                            "mount_path": os.path.join(self.target_outputs_dir, output_file.logical_name),
                            "host_path": output_file.host_path,
                            "read_only": False,
                        }
                     for output_file in request.output_files
                 ] + [
                        {
                            "name": "exec-script",
                            "mount_path": request.target_path,
                            "host_path": exec_script_path,
                            "read_only": True,
                        }
                 ]
        volumeMounts = [{"name": i['name'], "mount_path": i['mount_path'], "read_only": i['read_only']} for i in mounts]
        volumes = [{"name": i['name'], "host_path": {"path": i['host_path']}} for i in mounts]

        environment = request.environment + [
            f"INPUTS_DIR={self.target_inputs_dir}",
            f"OUTPUTS_DIR={self.target_outputs_dir}",
            f"ENFORCE_OUTPUTS_EXIST=1",
        ]

        script_path = request.target_path
        # volumes_from = []
        # if self._is_self_container():
        #     volumes_from += [ os.environ['HOSTNAME'] ]
        #     script_path = exec_script_path

        # Could set options like cpu_quote, mem_limit, restrict network access,
        # etc, to further lockdown the task.
        params = {
            "volumes": volumes,
        }

        container_params = {
            'env': [
                dict(zip( ("name","value"), i.split("=") ))
                for i in environment
            ],
            "volumneMounts": volumeMounts
        }

        return self.kubernetes_run(
            request.image, self._cmd(request, script_path),
            container_params,
            params
        )

    def reserve_exec_dir(self) -> ContextManager[str]:
        if not os.path.exists(self.exec_dir):
            os.mkdir(self.exec_dir)

        return tempfile.TemporaryDirectory(dir=self.exec_dir)

    def _cmd(self, request: RunRequest, exec_script_path: str):
        if request.run_with("node"):
            exec_dir = exec_script_path.replace(request.target_file, "")
            return ["/bin/bash", "-c", f"cd {exec_dir}; cp /app/package.json package.json; npm install; node {exec_script_path}"]
        elif request.run_with("r"):
            return ["Rscript", exec_script_path]
        else:
            return ["poetry", "run", "python", exec_script_path]

LocalProcess = Tuple[subprocess.Popen, IO[bin], str]


class LocalIsolator(Isolator[LocalProcess]):
    root_dir: str

    def is_running(self, t: LocalProcess) -> bool:
        p, io, f = t
        p.poll()
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
        inputs_dir = os.path.join(os.path.dirname(exec_script_path), "inputs")
        outputs_dir = os.path.join(os.path.dirname(exec_script_path), "outputs")
        stderr_log = os.path.join(os.path.dirname(exec_script_path), 'stderr')

        os.mkdir(inputs_dir)
        os.mkdir(outputs_dir)

        for input_file in request.input_files:
            in_path = os.path.abspath(input_file.host_path)
            isolated_path = os.path.join(inputs_dir, input_file.logical_name)
            print("Linking input: ", in_path, isolated_path, file=sys.stderr)

            os.symlink(in_path, isolated_path)

        for output_file in request.output_files:
            out_path = os.path.abspath(output_file.host_path)
            isolated_path = os.path.join(outputs_dir, output_file.logical_name)
            print("Linking output: ", out_path, isolated_path, file=sys.stderr)

            os.symlink(out_path, isolated_path)
            with open(out_path, "wb") as f:
                f.write(b'')

        print("Creating links for inputs and outputs", inputs_dir, outputs_dir, file=sys.stderr)

        environment = request.environment + [
            f"INPUTS_DIR={inputs_dir}",
            f"OUTPUTS_DIR={outputs_dir}",
            f"ENFORCE_OUTPUTS_EXIST=1",
        ]

        stderr_log_io = open(stderr_log, 'wb')
        return subprocess.Popen(
            [self._cmd(request), exec_script_path],
            stderr=stderr_log_io,
            env={k: v for k, v in [s.split('=', maxsplit=1) for s in environment]},
        ), stderr_log_io, stderr_log

    def _cmd(self, request: RunRequest):
        if request.run_with("python"):
            return shutil.which('python')
        else:
            # If we want to support alternate language  engines for local development,
            #   we'll have to pull the binary and archimedes-node
            #   into the archimedes image via a multi-stage docker build.
            # Otherwise, you just have to test via Vulcan.
            raise ValueError(
                "Script language is not supported in the local isolator. Use a docker isolator with an appropriate image."
            )


def run(request: RunRequest, isolator: Isolator[T], timeout = 60 * 60, remove = True) -> RunResult:
    res = RunResult(status='running', error=None)
    print(f"Preparing to run script using {request.isolator} isolator", file=sys.stderr)

    with isolator.reserve_exec_dir() as exec_dir:
        script_path = os.path.join(exec_dir, request.target_file)
        with open(script_path, "w") as script_file:
            script_file.write(request.script)

        print("Validating script...", file=sys.stderr)
        if request.run_with('python') and not run_checker([script_path]):
            raise ValueError(
                "Script did not conform to the archimedes dsl requirements."
            )

        print("Starting script...", file=sys.stderr)
        process: T = isolator.start(request, script_path)
        try:
            # Give cells up to 60 minutes.  In the future, with a proper async executor we would want to allow up to
            # a much larger amount of time.
            done = time.time() + timeout

            print("Waiting for process to finish...", file=sys.stderr)
            while isolator.is_running(process) and time.time() < done:
                time.sleep(1)

            timeout_reached = False
            if isolator.is_running(process):
                print("Timeout reached, forcing stop", file=sys.stderr)
                timeout_reached = True
                isolator.stop(process)

            code = isolator.wait(process)

            if isinstance(code, tuple):
                code, log = code
                res.status = 'failed'
                if code == 137:
                # KubernetesIsolator
                    if timeout_reached:
                        res.error = "Timeout reached, script stopped."
                    else:
                        res.error = "Out of memory."
                else:
                    # KubernetesIsolator outputs stdout and stderr directly
                    res.error = log
            elif code != 0:
                res.status = 'failed'
                if code == 137:
                    # When the container exits because of a 137,
                    #   stderr is an empty string, so not very
                    #   informative for debugging purposes.
                    if timeout_reached:
                        res.error = "Timeout reached, script stopped."
                    else:
                        res.error = "Out of memory."
                else:
                    res.error = isolator.get_stderr(process)
            else:
                res.status = 'done'
        finally:
            if remove and isinstance(process, Container):
                print("Cleaning up docker container", file=sys.stderr)
                process.remove()
            if remove and isinstance(process, str):
                print("Cleaning up kubernetes pod", file=sys.stderr)
                isolator.delete_pod(process, maybe_repeat=True)

    return res


def make_storage_file(s: str) -> StorageFile:
    parts = s.split(':', maxsplit=1)
    if len(parts) != 2:
        raise ValueError('files must be of the form name:/path/on/host')

    return StorageFile(
        logical_name=parts[0],
        host_path=parts[1],
    )


def main():
    import sys
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--file', help="The script file to run")
    parser.add_argument('--isolator', default='local', choices=['docker', 'local', 'kubernetes'])
    parser.add_argument('--image', default='etnaagent/archimedes:latest')
    parser.add_argument('--input', dest='inputs', default=[], action='append', help="input files of the form name:/path/on/host")
    parser.add_argument('--output', dest='outputs', default=[], action='append', help="output files of the form name:/path/on/host")
    parser.add_argument('-e', '--env', dest='env', action='append', help="environment variables of the form ABC=abc")
    parser.add_argument('--interpreter', default='python', choices=['python', 'node', 'r'])

    args = parser.parse_args()

    request: RunRequest = RunRequest(
        isolator=args.isolator,
        input_files=[make_storage_file(s) for s in args.inputs],
        output_files=[make_storage_file(s) for s in args.outputs],
        environment=args.env or [],
        script=(args.file and open(args.file, 'r').read()) or sys.stdin.read(),
        image=args.image,
        interpreter=args.interpreter
    )

    result = RunResult(status='done', error=f"Did not run, isolator {request.isolator} unrecognized")
    if request.isolator == 'docker':
        result = run(request, DockerIsolator())
    
    if request.isolator == 'kubernetes':
        result = run(request, KubernetesIsolator())

    if request.isolator == 'local':
        result = run(request, LocalIsolator())

    if result.error:
        print(result.error, file=sys.stderr)
    print(RunResult.schema().dumps(result))
