import hashlib
import logging
import queue
import time
from dataclasses import dataclass
from datetime import timedelta, datetime
from typing import Optional, List, Dict, Callable, Any, Iterator, Set, Mapping

import docker.errors
from airflow.exceptions import AirflowException
from airflow.models import TaskInstance, BaseOperator
from airflow.providers.docker.operators.docker import DockerOperator
from airflow.utils.strings import get_random_string
from dateutil import parser
from docker import types, APIClient
from docker.types import ConfigReference, Mount, SecretReference, Resources


@dataclass
class SwarmSharedData:
    data: bytes
    remote_path: str


def join_labels(labels: Dict[str, str]) -> List[str]:
    return [f"{k}={v}" for k, v in labels.items()]


@dataclass
class SwarmServiceDefinition:
    image: str
    command: List[str]
    mounts: Optional[List[Mount]]
    env: List[str]
    user: Optional[str]
    tty: bool
    configs: List[ConfigReference]
    secrets: List[SecretReference]
    resources: Optional[Resources]
    networks: List[types.NetworkAttachmentConfig]
    placement: types.Placement


def find_service(cli: APIClient, service_name: str) -> Dict:
    return cli.inspect_service(service_name, insert_defaults=True)


def find_local_network_ids(cli: APIClient, data: Dict) -> Set[str]:
    spec = data["Spec"]
    namespace = spec.get("Labels", {}).get("com.docker.stack.namespace", None)

    if namespace is not None:
        networks = cli.networks(
            filters=dict(label=join_labels({"com.docker.stack.namespace": namespace}))
        )
        return {n["Id"] for n in networks}

    return set()


def create_service_definition(
    data: Dict,
    local_network_ids: Optional[Set[str]] = None,
) -> SwarmServiceDefinition:

    spec = data["Spec"]
    task_template = spec["TaskTemplate"]
    container_spec = task_template["ContainerSpec"]

    networks = list(task_template.get("Networks", []))
    if local_network_ids is not None:
        networks = [n for n in networks if n["Target"] in local_network_ids]

    return SwarmServiceDefinition(
        image=container_spec["Image"],
        command=container_spec.get("Command", None),
        mounts=container_spec.get("Mounts", []),
        env=container_spec.get("Env", []),
        user=container_spec.get("User", None),
        tty=container_spec.get("TTY", None),
        configs=container_spec.get("Configs", []),
        secrets=container_spec.get("Secrets", []),
        resources=task_template.get("Resources", None),
        networks=networks,
        placement=task_template.get("Placement", None),
    )


def create_service_from_definition(
    cli: APIClient,
    service_definition: SwarmServiceDefinition,
    service_name: str,
    extr_env: Mapping[str, str],
    labeling: Dict[str, str],
):

    return cli.create_service(
        types.TaskTemplate(
            container_spec=types.ContainerSpec(
                image=service_definition.image,
                command=service_definition.command,
                mounts=service_definition.mounts,
                env=service_definition.env + list(f"{k}={v}" for k, v in extr_env.items()),
                user=service_definition.user,
                tty=service_definition.tty,
                configs=service_definition.configs,
                secrets=service_definition.secrets,
            ),
            restart_policy=types.RestartPolicy(condition="none"),
            resources=service_definition.resources,
            networks=service_definition.networks,
            placement=service_definition.placement,
        ),
        name=service_name,
        labels=labeling,
    )


def swarm_cleanup(cli: APIClient):
    for service in cli.services(
        filters={"label": join_labels(DockerSwarmOperator.service_labeling())}
    ):
        task_status = cli.tasks(filters={"service": service["ID"]})[0]["Status"]
        if task_status["State"] in [
            "complete",
            "failed",
            "shutdown",
            "rejected",
            "orphaned",
            "remove",
        ]:
            timestamp: str = task_status["Timestamp"]
            t: float = time.mktime(parser.parse(timestamp).timetuple())
            now = time.mktime(datetime.utcnow().timetuple())
            if t < now - timedelta(hours=2).total_seconds():
                cli.remove_service(service['ID'])

    for config in cli.configs(
            filters={"label": join_labels(DockerSwarmOperator.shared_data_labeling())}
    ):
        try:
            cli.remove_config(config['ID'])
        except docker.errors.APIError as e:
            # Config that happens to be shared by still running service will fail with 400,
            # but this ok. the cleanup will happen later.
            if e.response.status_code == 400:
                continue
            raise e

class DockerSwarmOperator(BaseOperator):
    swarm_shared_data: List[SwarmSharedData]
    command: List[str]
    serialize_last_output: Optional[Callable[[bytes], Any]] = None
    include_external_networks: Optional[bool]
    source_service: str
    cli: APIClient
    ti: TaskInstance
    terminated_service_state: Optional[str] = None
    docker_base_url: str
    env: Mapping[str, str]

    def __init__(
        self,
        *args,
        source_service: str,
        command: List[str],
        include_external_networks: Optional[bool] = False,
        swarm_shared_data: Optional[List[SwarmSharedData]] = None,
        serialize_last_output: Optional[Callable[[bytes], Any]] = None,
        docker_base_url="unix://var/run/docker.sock",
        env: Mapping[str, str] = dict(),
        **kwds,
    ):
        self.swarm_shared_data = swarm_shared_data or []
        self.serialize_last_output = serialize_last_output
        self.source_service = source_service
        self.include_external_networks = include_external_networks
        self.command = command
        self.docker_base_url = docker_base_url
        self.env = env
        super(DockerSwarmOperator, self).__init__(*args, **kwds)

    def execute(self, context) -> None:
        self.cli = APIClient(base_url=self.docker_base_url)
        self.ti = context["ti"]

        service_data = find_service(self.cli, self.source_service)

        local_network_ids = None
        if not self.include_external_networks:
            local_network_ids = find_local_network_ids(self.cli, service_data)

        command = DockerOperator.format_command(self.command)
        service_definition = create_service_definition(
            service_data,
            local_network_ids,
        )
        service_definition.command = command

        return self._run_service(service_definition)

    @staticmethod
    def shared_data_labeling():
        return {"etna.operators.swarm_operator.shared_data": "true"}

    @staticmethod
    def service_labeling():
        return {"etna.operators.swarm_operator.service": "true"}

    def _service_name(self) -> str:
        digest = hashlib.md5(f"{self.ti.dag_id}-{self.ti.task_id}-{self.ti.run_id}".encode('utf8')).hexdigest()
        return f"swarm_operator_{digest}"

    def _prepare_shared_data(self):
        configs: List[ConfigReference] = []

        for shared_data in self.swarm_shared_data:
            config = self.cli.create_config(
                f"{self.dag_id}-{self.task_id}-shared-data-{get_random_string()}",
                shared_data.data,
                self.shared_data_labeling(),
            )

            configs.append(
                ConfigReference(
                    config.id, config.name, filename=shared_data.remote_path
                )
            )

        return configs

    def _find_or_create_service(
        self, service_definition: SwarmServiceDefinition
    ) -> Dict:
        for service in self.cli.services(filters={"name": self._service_name()}):
            self.log.info("Attaching to running service: %s", str(self._service_name()))
            return service

        service = create_service_from_definition(
            self.cli, service_definition, self._service_name(), self.env, self.service_labeling()
        )

        self.log.info("Service started: %s", str(self._service_name()))
        return self.cli.inspect_service(service['ID'])

    def on_kill(self) -> None:
        if self.cli is not None:
            self.log.info("Removing docker service: %s", self.service["ID"])
            self.cli.remove_service(self.service["ID"])



    def _run_service(self, service_definition: SwarmServiceDefinition) -> Any:
        self.log.info("Starting docker service from image %s", service_definition.image)
        if not self.cli:
            raise Exception("The 'cli' should be initialized before!")

        self.service = self._find_or_create_service(service_definition)
        try:
            # wait for the service to start the task
            while not self.cli.tasks(filters={"service": self.service["ID"]}):
                self.log.info('Awaiting task to come up...')
                time.sleep(1)

            self.log.info('Consuming service logs now:')
            result = self._consume_logs()
            self.log.info("Service terminated: %s", self.terminated_service_state)

            if self.terminated_service_state != "complete":
                raise AirflowException(
                    "Service did not complete: " + repr(self.service)
                )
            return result
        finally:
            self._cleanup_service()

    def _cleanup_service(self):
        self.cli.remove_service(self.service["ID"])
        configs_json = self.service["Spec"]["TaskTemplate"]['ContainerSpec'].get("Configs", [])
        for config_spec_json in configs_json:
            config_json = self.cli.inspect_config(config_spec_json["ID"])
            config_labels = config_json["Spec"]["Labels"]

            if all(
                    config_labels.get(key) == v
                    for key, v in DockerSwarmOperator.shared_data_labeling().items()
            ):
                self.cli.remove_config(config_json["ID"])

    def _consume_logs(self) -> Any:
        logs_iter = consume_logs(self.cli, self.service['ID'], self.log)
        err = None

        # await the task being complete
        while True:
            # While consuming logs
            next(logs_iter)

            tasks = self.cli.tasks(filters={"service": self.service['ID']})
            status  = tasks[0]["Status"]
            err = status.get('Err')
            state = status["State"]
            if state in [
                "complete",
                "failed",
                "shutdown",
                "rejected",
                "orphaned",
                "remove",
            ]:
                self.terminated_service_state = state
                break
            time.sleep(1)

        # Then consume the remaining logs after a task completes, keeping the last log line
        last_line: bytes = next(logs_iter)

        if err:
            self.log.error(err)

        if self.terminated_service_state == 'complete' and self.serialize_last_output:
            return self.serialize_last_output(last_line)
        else:
            try:
                self.log.info(last_line.decode())
            except UnicodeDecodeError:
                self.log.info(last_line)

def await_status(
        cli: APIClient,
        service: Dict,
        result_queue: queue.Queue
):
    while True:
        state = cli.tasks(filters={"service": service['ID']})[0]["Status"]["State"]
        if state in [
            "complete",
            "failed",
            "shutdown",
            "rejected",
            "orphaned",
            "remove",
        ]:
            result_queue.put(state)
            break
        time.sleep(1)

# Addresses issue with current python client's service logs method.
# We force the entire logs to be consumed in a singular buffer rather than streamed.
# This isn't ideal but the streaming solution offered in the underlying implementation
# has many problems.  This should be ok since we still grab chunks using the 'since'
# property and are capable of consuming logs "as they are produced".  Worst case scenario,
# an extremely prolific logging service might cause network and memory bloating with
# the buffered approach, but this would be indicative of a deeper problem.
# In future maybe we can rewrite the python docker library's streaming implementation
# to fix underlying issues:
#   1. ability to cancel streaming (currently it removes socket timeout and does not allow retrieval
#      of underlying socket to cancel it, resulting in a leaking stream when services are completed)
#   2. it attempts to hack into internals of requests implementation to retrieve the underlying socket,
#      which is not portable when using the vcr requests library, nor is it portable when using an http
#      adapters or even different python versions.  This is a flaw in the way python http is implemented,
#      its attempt to hide underlying sockets causes grief (ie: stop hiding resources that are not actually private)
def hacky_service_logs(cli: APIClient, service_id: str, since: int, is_tty: bool) -> bytes:
    params = {
        'details': False,
        'follow': False,
        'stdout': True,
        'stderr': True,
        'since': since,
        'timestamps': True,
        'tail': 'all'
    }

    url = cli._url('/services/{0}/logs', service_id)
    res = cli._get(url, params=params, stream=False)
    return cli._get_result_tty(False, res, is_tty)

def consume_logs(
        cli: APIClient,
        service_id: str,
    log: logging.Logger,
) -> Iterator[bytes]:
    since = [0]
    last_line: List[bytes] = []

    is_tty = cli.inspect_service(
        service_id
    )['Spec']['TaskTemplate']['ContainerSpec'].get('TTY', False)

    def process_chunk(buff: bytes) -> bytes:
        lines = buff.split(b'\n')
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
        return b''

    while True:
        yield process_chunk(hacky_service_logs(cli, service_id, since[0], is_tty))