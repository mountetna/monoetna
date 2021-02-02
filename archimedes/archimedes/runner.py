import os.path
import os
import threading
from dataclasses import dataclass
from hashlib import md5
from typing import List, Iterable
from docker import DockerClient
from docker.types import Mount
from archimedes.checker import run as run_checker
import tempfile


@dataclass
class RunResult:
    status: str
    error: List[str]

class Runner:
    def open_key(self, key, mode):
        # TODO, this is just a demonstration.  Ideally this would be a input/output stream to
        # external storage service.
        if mode.startswith("r"):
            return open(os.path.join(os.path.dirname(__file__), "..", "tests", "fixtures", "runner_inputs", key), mode)
        return tempfile.TemporaryFile(mode)

    # A simple entrypoint demonstrating a process to run.
    def run(self, script_text, image="archimedes-base", output_keys: Iterable[bin] = tuple(), input_keys: Iterable[bin] = tuple(), environment: Iterable[str] = tuple()):
        """
        :param script_text: The python script's text contents to run.
        :param output_keys: Keys to store process outputs. Any iterable.
        :param input_keys:  Keys to find process inputs. Any iterable.
        :param environment: Iterable of strings in the form "MYENVKEY=myvalue"
        :return: For now, nothing, but in the future, this function may become asynchronous and
            attach to a scheduler, in which case it would return the status of the job being scheduled.
        """
        with tempfile.TemporaryDirectory() as inputs_dir, tempfile.TemporaryDirectory() as outputs_dir, tempfile.TemporaryDirectory() as exec_dir:
            script_path = os.path.join(exec_dir, "script.py")
            with open(script_path, "w") as script_file:
                script_file.write(script_text)

            if not run_checker([script_path]):
                raise ValueError("Script did not conform to the archimedes dsl requirements.")

            mounts = [
                Mount(
                    target="/var/run/vulcan/exec",
                    source=exec_dir,
                    read_only=True
                ),
                Mount(
                    target="/var/run/vulcan/inputs",
                    source=inputs_dir,
                    read_only=True
                ),
                Mount(
                    target="/var/run/vulcan/outputs",
                    source=outputs_dir,
                    read_only=False
                ),
            ]

            output_pipes = [
                (os.mkfifo(os.path.join(outputs_dir, md5(output_key).hexdigest())), output_key)
                for output_key in output_keys
            ]

            input_pipes = [
                (os.mkfifo(os.path.join(outputs_dir, md5(input_key).hexdigest())), input_key)
                for input_key in output_keys
            ]

            def consume_outputs():
                pass

            def produce_inputs():
                pass

            consume_outputs_thread = threading.Thread(target=consume_outputs)
            consume_outputs_thread.start()

            produce_inputs_thread = threading.Thread(target=produce_inputs)
            produce_inputs_thread.start()

            docker_cli: DockerClient = DockerClient.from_env()
            # Could set options like cpu_quote, mem_limit, restrict network access,
            # etc, to further lockdown the task.
            docker_cli.containers.run(image,
                                      "poetry run /var/run/vulcan/exec/script.py",
                                      remove=True,
                                      environment=list(environment),
                                      # Arbitrary limit of 2G for entire image size.
                                      storage_opt=dict(size="2G"),
                                      mounts=mounts)
