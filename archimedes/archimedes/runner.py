import os.path
import os
import shutil
from concurrent.futures import Future
from dataclasses import dataclass
from typing import List, Iterable, Tuple, Dict, Any, BinaryIO, Set
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED

import requests
from docker import DockerClient
from docker.types import Mount
from docker.models.containers import Container
from archimedes.checker import run as run_checker
import tempfile


@dataclass
class RunResult:
    status: str
    error: List[str]


class Runner:
    def open_key(self, key, mode) -> BinaryIO:
        # TODO, this is just a demonstration.  Ideally this would be a input/output stream to
        # external storage service.
        if mode.startswith("r"):
            return open(
                os.path.join(
                    os.path.dirname(__file__),
                    "..",
                    "tests",
                    "fixtures",
                    "runner_inputs",
                    key,
                ),
                "rb",
            )
        # Trick the type checker.
        _result: Any = tempfile.TemporaryFile("wb")
        return _result

    def _prepare_fifo_files_from_keys(
        self, names_and_keys: Iterable[Tuple[bin, bin]], keys_dir
    ):
        for name, key in names_and_keys:
            path = os.path.join(keys_dir, name)
            os.mkfifo(path)
            yield path, key

    def _prepare_container(
        self, image, inputs_dir, outputs_dir, exec_dir, environment: Iterable[str]
    ) -> Container:

        target_inputs_dir = "/var/run/vulcan/inputs"
        target_outputs_dir = "/var/run/vulcan/outputs"

        mounts = [
            Mount(target="/var/run/vulcan/exec", source=exec_dir, read_only=True),
            Mount(target=target_inputs_dir, source=inputs_dir, read_only=True),
            Mount(target=target_outputs_dir, source=outputs_dir, read_only=False),
        ]
        environment += (
            f"INPUTS_DIR={target_inputs_dir}",
            f"OUTPUTS_DIR={target_outputs_dir}",
            f"ENFORCE_OUTPUTS_EXIST=1",
        )

        docker_cli: DockerClient = DockerClient.from_env()
        # Could set options like cpu_quote, mem_limit, restrict network access,
        # etc, to further lockdown the task.
        return docker_cli.containers.run(
            image,
            "poetry run /var/run/vulcan/exec/script.py",
            auto_remove=True,
            detach=True,
            environment=list(environment),
            # Arbitrary limit of 2G for entire image size.
            storage_opt=dict(size="2G"),
            mounts=mounts,
        )

    # A simple entrypoint demonstrating a process to run.
    # TODO: When detaching containers, label them in such a way that they are reaped
    def run(
        self,
        script_text: str,
        image="vulcan_base",
        output_names_and_keys: Iterable[Tuple[bin, bin]] = tuple(),
        input_names_and_keys: Iterable[Tuple[bin, bin]] = tuple(),
        environment: Iterable[str] = tuple(),
    ):
        """
        :param script_text: The python script's text contents to run.
        :param output_names_and_keys: Iterable of cell input names and key paths
        :param input_names_and_keys:  Iterable of cell output names and key paths
        :param environment: Iterable of strings in the form "MYENVKEY=myvalue"
        :return: For now, nothing, but in the future, this function may become asynchronous and
            attach to a scheduler, in which case it would return the status of the job being scheduled.
        """
        with tempfile.TemporaryDirectory() as inputs_dir, tempfile.TemporaryDirectory() as outputs_dir, tempfile.TemporaryDirectory() as exec_dir:
            script_path = os.path.join(exec_dir, "script.py")
            with open(script_path, "w") as script_file:
                script_file.write(script_text)

            if not run_checker([script_path]):
                raise ValueError(
                    "Script did not conform to the archimedes dsl requirements."
                )

            output_pipes = list(
                self._prepare_fifo_files_from_keys(output_names_and_keys, outputs_dir)
            )
            input_pipes = list(
                self._prepare_fifo_files_from_keys(input_names_and_keys, inputs_dir)
            )

            container = self._prepare_container(
                image, inputs_dir, outputs_dir, exec_dir, environment
            )
            result = self._execute_container(container, output_pipes, input_pipes)
            if result["StatusCode"] != 1:
                print(result)
                raise RuntimeError("Script failed, check the logs")

    def _consume_output_key(self, src, key):
        with self.open_key(key, 'wb') as key_file:
            with open(src, "rb") as output_file:
                shutil.copyfileobj(output_file, key_file)

    def _produce_input_key(self, dest, key):
        with self.open_key(key, 'rb') as key_file:
            with open(dest, "wb") as input_file:
                shutil.copyfileobj(key_file, input_file)

    def _execute_container(
        self,
        container: Container,
        output_pipes: List[Tuple[str, str]],
        input_pipes: List[Tuple[str, str]],
    ) -> Dict[str, Any]:
        with ThreadPoolExecutor(
            max_workers=len(output_pipes) + len(input_pipes) + 1
        ) as executor:
            consuming: Set[Future[Any]] = set(
                executor.submit(self._consume_output_key, path, key)
                for path, key in output_pipes
            )
            producing: Set[Future[Any]] = set(
                executor.submit(self._produce_input_key, path, key)
                for path, key in input_pipes
            )

            def run_container():
                while True:
                    try:
                        return container.wait(timeout=1)
                    except requests.exceptions.ReadTimeout:
                        pass

            ran = executor.submit(run_container)

            try:
                done: Set[Future[Any]] = set()
                pending: Set[Future[Any]] = consuming.union(producing).union({ran})
                f: Future

                while ran not in done:
                    next_done, pending = wait(pending, return_when=FIRST_COMPLETED)
                    done.update(next_done)

                    for f in next_done:
                        if f is ran or f.exception():
                            break

                if ran.done():
                    pass

                for f in pending:
                    f.cancel()

                executor.shutdown()

                if ran.cancelled():
                    raise RuntimeError("Script failed, pipe was broken")

                if ran.exception():
                    raise ran.exception()
                return ran.result()
            finally:
                container.stop()
