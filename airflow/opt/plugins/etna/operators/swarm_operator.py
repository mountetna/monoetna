import time
from dataclasses import dataclass
from datetime import timedelta, datetime
from dateutil import parser
from typing import Optional, List, Dict, Callable, Any, Union, Iterator, Set

import docker.errors
from airflow.exceptions import AirflowException
from airflow.models import TaskInstance, BaseOperator
from airflow.providers.docker.operators.docker import DockerOperator
from airflow.utils.strings import get_random_string
from docker import types, APIClient

from airflow.providers.docker.operators.docker_swarm import (
    DockerSwarmOperator as OrigDockerSwarmOperator,
)
from docker.models.services import Service
from docker.types import ConfigReference, Mount, SecretReference, ServiceMode, Resources


@dataclass
class SwarmSharedData:
    data: bytes
    remote_path: str


def join_labels(labels: Dict[str, str]) -> List[str]:
    return [f"{k}=={v}" for k, v in labels.items()]


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
    task_template = spec["TaskTemplate"]
    namespace = task_template.get("Labels", {}).get("com.docker.stack.namespace", None)

    if namespace is not None:
        networks = cli.networks(
            filters=dict(labels=join_labels({"com.docker.stack.namespace": namespace}))
        )
        return {n["ID"] for n in networks}

    return set()


def create_service_definition(
    data: Dict,
    local_network_ids: Optional[Set[str]] = None,
) -> SwarmServiceDefinition:

    spec = data["Spec"]
    task_template = spec["TaskTemplate"]
    container_spec = task_template["ContainerSpec"]

    networks = task_template.get("Networks", [])
    if local_network_ids:
        networks = [n for n in networks if n["ID"] in local_network_ids]

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
    labeling: Dict[str, str],
):
    return cli.create_service(
        types.TaskTemplate(
            container_spec=types.ContainerSpec(
                image=service_definition.image,
                command=service_definition.command,
                mounts=service_definition.mounts,
                env=service_definition.env,
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


# TODO: Executor can hook this up.
def swarm_cleanup(cli: APIClient):
    for config in cli.configs(
        filters={"label": join_labels(DockerSwarmOperator.shared_data_labeling())}
    ):
        try:
            config.remove()
        except docker.errors.APIError:
            pass  # ignore, maybe being used still.

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
            t: datetime = parser.parse(timestamp)
            if t < datetime.now() - timedelta(hours=1):
                try:
                    service.remove()
                except docker.errors.APIError:
                    pass  # ignore, maybe being used still.


class DockerSwarmOperator(BaseOperator):
    swarm_shared_data: List[SwarmSharedData]
    command: List[str]
    serialize_last_output: Optional[Callable[[bytes], Any]] = None
    include_external_networks: Optional[bool]
    source_service: str
    cli: APIClient
    ti: TaskInstance

    def __init__(
        self,
        source_service: str,
        command: List[str],
        *args,
        include_external_networks: Optional[bool] = False,
        swarm_shared_data: Optional[List[SwarmSharedData]] = None,
        serialize_last_output: Optional[Callable[[bytes], Any]] = None,
        **kwds,
    ):
        self.swarm_shared_data = swarm_shared_data or []
        self.serialize_last_output = serialize_last_output
        self.source_service = source_service
        self.include_external_networks = include_external_networks
        self.command = command
        super(DockerSwarmOperator, self).__init__(*args, **kwds)

    def execute(self, context) -> None:
        self.cli = APIClient(base_url="unix://var/run/docker.sock")
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

    def _service_name(self):
        return f"{self.ti.dag_id}-{self.ti.task_id}-{self.ti.run_id}"

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
    ) -> Service:
        for service in self.cli.services(filters={"name": self._service_name()}):
            self.log.info("Attaching to running service: %s", str(self._service_name()))
            return service

        service = create_service_from_definition(
            self.cli, service_definition, self._service_name(), self.service_labeling()
        )

        self.log.info("Service started: %s", str(self._service_name()))
        return service

    def on_kill(self) -> None:
        if self.cli is not None:
            self.log.info("Removing docker service: %s", self.service["ID"])
            self.cli.remove_service(self.service["ID"])

    def _run_service(self, service_definition: SwarmServiceDefinition) -> None:
        self.log.info("Starting docker service from image %s", service_definition.image)
        if not self.cli:
            raise Exception("The 'cli' should be initialized before!")

        self.service = self._find_or_create_service(service_definition)

        # wait for the service to start the task
        while not self.cli.tasks(filters={"service": self.service["ID"]}):
            time.sleep(1)

        result = self._consume_logs(service_definition)

        while True:
            if self._has_service_terminated():
                self.log.info(
                    "Service status before exiting: %s", self._service_status()
                )
                break

        try:
            if self._service_status() != "complete":
                raise AirflowException(
                    "Service did not complete: " + repr(self.service)
                )
            return result
        finally:
            self._cleanup_service()

    def _cleanup_service(self):
        self.cli.remove_service(self.service["ID"])
        configs_json = self.service["Spec"]["TaskTemplate"]["Configs"]
        for config_spec_json in configs_json:
            config_json = self.cli.inspect_config(config_spec_json["ID"])
            config_labels = config_json["Spec"]["Labels"]

            if all(
                config_labels.get(key) == v
                for key, v in self.shared_data_labeling().items()
            ):
                self.cli.remove_config(config_json["ID"])

    def _service_status(self) -> Optional[str]:
        if not self.cli:
            raise Exception("The 'cli' should be initialized before!")
        return self.cli.tasks(filters={"service": self.service["ID"]})[0]["Status"][
            "State"
        ]

    def _has_service_terminated(self) -> bool:
        status = self._service_status()
        return status in [
            "complete",
            "failed",
            "shutdown",
            "rejected",
            "orphaned",
            "remove",
        ]

    def _consume_logs(self, service_definition: SwarmServiceDefinition) -> Any:
        log_chunk_iter = self.cli.service_logs(
            self.service["ID"],
            follow=True,
            stdout=True,
            stderr=True,
            is_tty=service_definition.tty,
        )
        last_line_buffer: List[bytes] = []

        for line in consume_lines(log_chunk_iter, last_line_buffer):
            self.log.info(line)

        if self.serialize_last_output:
            return self.serialize_last_output(b"".join(last_line_buffer))
        else:
            try:
                self.log.info(b"".join(last_line_buffer))
            except UnicodeDecodeError:
                pass


def consume_lines(
    chunk_iter: Iterator[bytes], last_line_buffer: List[bytes]
) -> Iterator[str]:
    buffer_ready = False

    for chunk in chunk_iter:
        if chunk == b"\n":
            buffer_ready = True
        else:
            if buffer_ready:
                a = b"".join(last_line_buffer)
                last_line_buffer.clear()
                buffer_ready = False

                try:
                    b = a.decode()
                    yield b
                except UnicodeDecodeError:
                    pass

            last_line_buffer.append(chunk)
