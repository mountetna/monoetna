import functools
import os.path
import json
import logging
import time
from dataclasses import dataclass
from typing import (
    List,
    Optional,
    Callable,
    Any,
    Mapping,
    Iterator,
    Set,
    Tuple,
    Union,
    MutableMapping,
)

import hashlib

import io

import tarfile
from airflow import AirflowException
from airflow.models import TaskInstance, XCom, BaseOperator
from airflow.models.xcom_arg import XComArg
from dateutil import parser
from docker import APIClient

from etna.utils.gater import RetryGater


@dataclass
class SwarmSharedData:
    data: bytes
    remote_path: str


class DockerOperatorBase(BaseOperator):
    terminated_state: Optional[str]
    _swarm_shared_data: List[SwarmSharedData]
    command: List[str]
    serialize_last_output: Optional[Callable[[bytes], Any]] = None
    cli: APIClient
    ti: TaskInstance
    docker_base_url: str
    env: Mapping[str, str]
    args: List[str]

    file_inputs: MutableMapping[str, Union[XComArg, str, bytes]]
    resolved_swarm_shared_data: List[SwarmSharedData]

    template_fields = ("env",)
    template_fields_renderers = {"env": "json"}

    def __init__(
        self,
        *args,
        command: List[str],
        swarm_shared_data: Optional[List[SwarmSharedData]] = None,
        serialize_last_output: Optional[Callable[[bytes], Any]] = None,
        docker_base_url="unix://var/run/docker.sock",
        env: Mapping[str, str] = dict(),
        **kwds,
    ):
        self._swarm_shared_data = swarm_shared_data or []
        self.serialize_last_output = serialize_last_output
        self.command = command
        self.docker_base_url = docker_base_url
        self.env = env
        self.file_inputs = {}
        self.resolved_swarm_shared_data = []
        self.args = []

        super(DockerOperatorBase, self).__init__(*args, **kwds)

    def accepts(self, *file_names: str):
        self.args = list(file_names)
        return self

    def __call__(self, *args):
        if len(args) != len(self.args):
            raise TypeError(
                f"docker operator should use .accepts to match inputs size, configured for {self.args} but got {len(args)}"
            )
        for xarg, file_arg in zip(args, self.args):
            self[file_arg] = xarg
        return XComArg(self)

    def pre_execute(self, context) -> None:
        super(DockerOperatorBase, self).pre_execute(context)
        self.cli = APIClient(base_url=self.docker_base_url)
        self.ti = context["ti"]

        self.resolved_swarm_shared_data = []
        for path, file_input in self.file_inputs.items():
            bin: bytes
            if isinstance(file_input, XComArg):
                bin = XCom.serialize_value(file_input.resolve(context))
            elif isinstance(file_input, bytes):
                bin = file_input
            elif isinstance(file_input, str):
                if not os.path.exists(file_input):
                    raise AirflowException(
                        f"docker operator input {path} was not a valid path: {file_input}."
                    )

                if os.path.isdir(file_input):
                    buff = io.BytesIO(b"")
                    tar = tarfile.open(fileobj=buff)
                    tar.add(file_input, arcname=os.path.basename(path))
                    tar.close()
                    bin = buff.read()
                else:
                    bin = open(file_input, "rb").read()
            else:
                raise f"Unexpected input to docker operator: Expected XComArg, str, or bytes, found {file_input.__class__.__name__}"

            self.resolved_swarm_shared_data.append(
                SwarmSharedData(remote_path=path, data=bin)
            )

    @property
    def swarm_shared_data(self):
        return self._swarm_shared_data + self.resolved_swarm_shared_data

    def set_input(self, key: str, value: XComArg) -> "DockerOperatorBase":
        self[key] = value
        return self

    def __setitem__(self, key, value):
        if (
            not isinstance(value, XComArg)
            and not isinstance(value, str)
            and not isinstance(value, bytes)
        ):
            raise AirflowException(
                f"'{key}' was not instance of XComArg, str, or bytes!"
            )
        if not isinstance(key, str):
            raise AirflowException(
                f"{self.__class__.__name__} only accepts string index keys"
            )
        if isinstance(value, XComArg):
            self.set_upstream(value)
        if not key.startswith("/") or len(key) < 2:
            raise AirflowException(
                f"Key {key} is an invalid docker operator argument, must be a absolute path."
            )
        self.file_inputs[key] = value

    @staticmethod
    def shared_data_labeling():
        return {"etna.operators.swarm_operator.shared_data": "true"}

    @property
    def task_name(self) -> str:
        digest = hashlib.md5(
            f"{self.ti.dag_id}-{self.ti.task_id}-{self.ti.run_id}".encode("utf8")
        ).hexdigest()
        return f"{self.__class__.__name__}_{digest}"

    @property
    def successful_states(self) -> Set[str]:
        raise NotImplementedError("Oh no!")

    @property
    def completed_states(self) -> Set[str]:
        raise NotImplementedError("Oh no!")

    @property
    def _next_log_batch_since(self) -> Callable[[int], bytes]:
        raise NotImplementedError("Oh no!")

    @property
    def cleanup(self):
        raise NotImplementedError("oh no!")

    def on_kill(self) -> None:
        self.cleanup()
        return super(DockerOperatorBase, self).on_kill()

    def _check_task(self) -> Tuple[str, Optional[str]]:
        raise NotImplementedError("oh no!")

    def _start_task(self):
        raise NotImplementedError("oh no!")

    def execute(self, context: Any):
        self._start_task()
        try:
            self.log.info("Consuming logs now:")
            result = self._process_task_run()
            self.log.info("Task terminated: %s", self.terminated_state)

            if not self.terminated_state in self.successful_states:
                raise AirflowException("Task did not complete successfully")
            return result
        finally:
            self.cleanup()

    def _process_task_run(self) -> Any:
        logs_iter = write_logs_and_yield_last(self._next_log_batch_since, self.log)
        err = None

        # await the task being complete
        while True:
            # While consuming logs
            next(logs_iter)

            state, err = self._check_task()
            if state in self.completed_states:
                self.terminated_state = state
                break
            time.sleep(1)

        # Then consume the remaining logs after a task completes, keeping the last log line
        last_line: bytes = next(logs_iter)

        if err:
            self.log.error(err)

        if (
            self.terminated_state in self.successful_states
            and self.serialize_last_output
        ):
            return self.serialize_last_output(last_line)
        else:
            try:
                self.log.info(last_line.decode())
            except UnicodeDecodeError:
                self.log.info(last_line)


LOG_LINE_BATCH_SIZE = 16384  # 16K
MAX_BATCH_SIZE = 500  # 16K * 500 ~= 8mb
LOG_TIMESTAMP_LEN = 30


def write_logs_and_yield_last(
    next_batch_since: Callable[[int], bytes],
    log: logging.Logger,
) -> Iterator[bytes]:
    since_s = b""
    since_t = 0
    last_line: Optional[bytes] = None
    next_batch_buff: List[bytes] = []

    gater = RetryGater(5, exp=1.5)

    while True:
        buff = next_batch_since(since_t)

        lines = buff.split(b"\n")
        rpc_error = False

        for line in lines:
            if not line:
                continue

            if line.startswith(b"Error"):
                # rpc error incoming.  Attempt a fresh new fetch. do NOT allow a yield after an error
                # since we need to ensure some sort of conclusion to fetching logs (and if we cannot, we'd rather
                # timeout or error from connection, not rpc issues).
                rpc_error = True
                gater.gate(
                    AirflowException(
                        f"Could not read docker logs, received rpc error: {line.decode('utf8')}"
                    )
                )
                break

            # We need a 'timestamp' for the since call, but the string associated with the time may have higher precision
            # that we care about.  We compare using the higher precision string, but pass to the api the int.  This
            # will result in some overlap of results but we can dedup using the string for comparison.
            # In theory, logs can be lost if even the string precision is not string enough, but there isn't an easy way
            # around this since the api lacks a strictly increasing non duplicated value for consuming logs.
            date_s = line[:LOG_TIMESTAMP_LEN]
            date = parser.parse(date_s)
            t = time.mktime(date.timetuple())

            if date_s <= since_s:
                continue

            since_t = t
            since_s = date_s
            line = line[31:]

            # log messages > 16K are broken into batches without newline separation, but including the timestamp
            # header which must be removed.
            if len(line) > LOG_LINE_BATCH_SIZE:
                line_segments = []
                offset = 0
                for i in range(0, MAX_BATCH_SIZE):
                    if offset > len(line):
                        break
                    line_segments.append(line[offset : offset + LOG_LINE_BATCH_SIZE])
                    offset += LOG_LINE_BATCH_SIZE + LOG_TIMESTAMP_LEN + 1
                line = b"".join(line_segments)

            next_batch_buff.append(line)

        if rpc_error:
            continue

        if next_batch_buff:
            if last_line is not None:
                try:
                    log.info(last_line.decode())
                except UnicodeDecodeError:
                    log.info(last_line)

            for line in next_batch_buff[:-1]:
                try:
                    log.info(line.decode())
                except UnicodeDecodeError:
                    log.info(line)

            last_line = next_batch_buff[-1]
            next_batch_buff = []

        gater.reset()
        if last_line is not None:
            yield last_line
        else:
            yield b""
