import logging
import time
from typing import List, Optional, Callable, Any, Mapping, Iterator, Set, Tuple

import hashlib

from airflow import AirflowException
from airflow.models import BaseOperator, TaskInstance
from dateutil import parser
from docker import APIClient

from etna import SwarmSharedData


class DockerOperatorBase(BaseOperator):
    terminated_state: Optional[str]
    swarm_shared_data: List[SwarmSharedData]
    command: List[str]
    serialize_last_output: Optional[Callable[[bytes], Any]] = None
    cli: APIClient
    ti: TaskInstance
    docker_base_url: str
    env: Mapping[str, str]

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
        self.swarm_shared_data = swarm_shared_data or []
        self.serialize_last_output = serialize_last_output
        self.command = command
        self.docker_base_url = docker_base_url
        self.env = env
        super(DockerOperatorBase, self).__init__(*args, **kwds)

    def pre_execute(self, context) -> None:
        super(DockerOperatorBase, self).pre_execute(context)
        self.cli = APIClient(base_url=self.docker_base_url)
        self.ti = context["ti"]

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


def write_logs_and_yield_last(
    next_batch_since: Callable[[int], bytes],
    log: logging.Logger,
) -> Iterator[bytes]:
    since = [0]
    last_line: List[bytes] = []

    def process_chunk(buff: bytes) -> bytes:
        lines = buff.split(b"\n")
        results: List[bytes] = []
        for line in lines:
            if not line:
                continue
            date = parser.parse(line[:30])
            t = time.mktime(date.timetuple())

            # There isn't a way to consume 'since' with a exclusive >, meaning that
            # between calls to the logs we get overlapping lines on the last since
            # value passed.  We compare the time and drop those.
            if t <= since[0]:
                continue

            since[0] = t
            line = line[31:]
            results.append(line)

        if results:
            if last_line:
                try:
                    log.info(last_line[0].decode())
                except UnicodeDecodeError:
                    log.info(last_line[0])
            last_line.clear()

            for line in results[:-1]:
                try:
                    log.info(line.decode())
                except UnicodeDecodeError:
                    log.info(line)

            last_line.append(results[-1])

        if last_line:
            return last_line[0]
        return b""

    while True:
        yield process_chunk(next_batch_since(since[0]))
