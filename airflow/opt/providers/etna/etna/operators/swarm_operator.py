import hashlib
import time
from dataclasses import dataclass
from datetime import timedelta, datetime
from typing import Optional, List, Dict, Callable, Set, Mapping, Tuple

import requests.exceptions
import docker.errors
from airflow import AirflowException
from airflow.providers.docker.operators.docker import DockerOperator
from airflow.utils.strings import get_random_string
from dateutil import parser
from docker import types, APIClient
from docker.types import ConfigReference, Mount, SecretReference, Resources

from etna.operators.docker_operator_base import DockerOperatorBase
from etna.utils.gater import RetryGater


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
    allow_manager_nodes: Optional[bool] = False,
    resources: Optional[Resources] = None,
    additional_constraints: Optional[List[str]] = None,
) -> SwarmServiceDefinition:

    spec = data["Spec"]
    task_template = spec["TaskTemplate"]
    container_spec = task_template["ContainerSpec"]

    networks = list(task_template.get("Networks", []))
    if local_network_ids is not None:
        networks = [n for n in networks if n["Target"] in local_network_ids]

    placement = task_template.setdefault("Placement", {})
    if not allow_manager_nodes:
        constraints = [
            c for c in placement.get("Constraints", [])
            if "node.role" not in c
        ]
        constraints.append("node.role!=manager")
        placement["Constraints"] = constraints

    if additional_constraints is not None:
        placement["Constraints"] += additional_constraints

    return SwarmServiceDefinition(
        image=container_spec["Image"],
        command=container_spec.get("Command", None),
        mounts=container_spec.get("Mounts", []),
        env=container_spec.get("Env", []),
        user=container_spec.get("User", None),
        tty=container_spec.get("TTY", None),
        configs=container_spec.get("Configs", []),
        secrets=container_spec.get("Secrets", []),
        resources=resources or container_spec.get("Resources", None),
        networks=networks,
        placement=placement,
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
                env=service_definition.env
                + list(f"{k}={v}" for k, v in extr_env.items()),
                user=service_definition.user,
                tty=service_definition.tty,
                configs=service_definition.configs,
                secrets=service_definition.secrets,
            ),
            log_driver=types.DriverConfig(
                "local"
            ),  # Ensures logs are easily recoverable.
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
                cli.remove_service(service["ID"])

    for config in cli.configs(
        filters={"label": join_labels(DockerSwarmOperator.shared_data_labeling())}
    ):
        try:
            cli.remove_config(config["ID"])
        except docker.errors.APIError as e:
            # Config that happens to be shared by still running service will fail with 400,
            # but this ok. the cleanup will happen later.
            if e.response.status_code == 400:
                continue
            raise e


class DockerSwarmOperator(DockerOperatorBase):
    include_external_networks: Optional[bool]
    allow_manager_nodes: Optional[bool]
    source_service: str
    service: Mapping
    # NanoCpus (1 cpu = 1e9)
    cpu_limit: Optional[int]
    # Bytes (1KB = 1024, 1MB = 1024 * 1024, 1GIG = 1024 * 1024 * 1024)
    mem_limit: Optional[int]
    cpu_reservation: Optional[int]
    mem_reservation: Optional[int]
    additional_constraints: Optional[List[str]]

    def __init__(
            self,
        *args,
        source_service: str,
        include_external_networks: Optional[bool] = False,
        allow_manager_nodes: Optional[bool] = False,
        # NanoCpus (1 cpu = 1e9)
        cpu_limit: Optional[int] = int(1e9 * 0.5),
        # Bytes (1KB = 1024, 1MB = 1024 * 1024, 1GIG = 1024 * 1024 * 1024)
        mem_limit: Optional[int] = None,
        cpu_reservation: Optional[int] = None,
        mem_reservation: Optional[int] = None,
        additional_constraints: Optional[List[str]] = None,
        **kwds,
    ):
        self.mem_reservation = mem_reservation
        self.cpu_reservation = cpu_reservation
        self.mem_limit = mem_limit
        self.cpu_limit = cpu_limit
        self.source_service = source_service
        self.include_external_networks = include_external_networks
        self.allow_manager_nodes = allow_manager_nodes
        self.additional_constraints = additional_constraints
        super(DockerSwarmOperator, self).__init__(*args, **kwds)

    def _start_task(self):
        service_data = find_service(self.cli, self.source_service)
        local_network_ids = None
        if not self.include_external_networks:
            local_network_ids = find_local_network_ids(self.cli, service_data)

        command = DockerOperator.format_command(self.command)
        service_definition = create_service_definition(
            service_data,
            local_network_ids,
            allow_manager_nodes=self.allow_manager_nodes,
            additional_constraints=self.additional_constraints,
            resources=Resources(
                cpu_limit=self.cpu_limit,
                cpu_reservation=self.cpu_reservation,
                mem_limit=self.mem_limit,
                mem_reservation=self.mem_reservation,
            )
        )
        service_definition.command = command

        for config in self._prepare_shared_data():
            service_definition.configs.append(config)

        self.service = self._find_or_create_service(service_definition)

        while not self._get_task():
            self.log.info("Awaiting task to come up...")
            time.sleep(1)

    def _get_task(self) -> Optional[Mapping]:
        tasks = self.cli.tasks(filters={"service": self.service["ID"]})
        if tasks:
            return tasks[0]
        return None

    def _check_task(self) -> Tuple[str, Optional[str]]:
        task = self._get_task()
        if task:
            status = task["Status"]
            err = status.get("Err")
            state = status["State"]
            return state, err
        return "starting", None

    @property
    def _next_log_batch_since(self) -> Callable[[int], bytes]:
        is_tty = self.cli.inspect_service(self.service["ID"])["Spec"]["TaskTemplate"][
            "ContainerSpec"
        ].get("TTY", False)
        task = self._get_task()

        def next_log_batch(since: int):
            if task is None:
                # _start_task should wait until the task exists before proceeding to invoke the next_log_batch_since.
                # We also do NOT refresh the task after the initial fetch -- the task could legitimately go away while
                # reading, but we only care about the historical task identifier.
                raise AirflowException(
                    "Could not logs of service, task was not available before reading logs."
                )

            gater = RetryGater(5)

            while True:
                try:
                    return hacky_task_logs(self.cli, task["ID"], since, is_tty)
                except requests.exceptions.HTTPError as e:
                    if e.response.status_code >= 500:
                        gater.gate(e)
                        continue
                    raise e
                except requests.exceptions.ConnectionError as e:
                    gater.gate(e)
                    continue
                except requests.exceptions.Timeout as e:
                    gater.gate(e)
                    continue

        return next_log_batch

    @property
    def successful_states(self) -> Set[str]:
        return {
            "complete",
        }

    @property
    def completed_states(self) -> Set[str]:
        return self.successful_states.union(
            {
                "failed",
                "shutdown",
                "rejected",
                "orphaned",
                "remove",
            }
        )

    @staticmethod
    def shared_data_labeling():
        return {"etna.operators.swarm_operator.shared_data": "true"}

    @staticmethod
    def service_labeling():
        return {"etna.operators.swarm_operator.service": "true"}

    def _prepare_shared_data(self):
        configs: List[ConfigReference] = []

        for shared_data in self.swarm_shared_data:
            config_name = (
                f"{self.dag_id}-{self.task_id}-shared-data-{get_random_string()}"
            )

            if len(config_name) > 64:
                config_name = hashlib.md5(config_name.encode('utf-8')).hexdigest()

            config = self.cli.create_config(
                config_name,
                shared_data.data,
                self.shared_data_labeling(),
            )

            configs.append(
                ConfigReference(
                    config["ID"], config_name, filename=shared_data.remote_path
                )
            )

        return configs

    def _find_or_create_service(
        self, service_definition: SwarmServiceDefinition
    ) -> Dict:
        for service in self.cli.services(filters={"name": self.task_name}):
            self.log.info("Attaching to running service: %s", str(self.task_name))
            return service

        self.log.info(
            "Starting %s inside of service %s", self.command, self.source_service
        )
        service = create_service_from_definition(
            self.cli,
            service_definition,
            self.task_name,
            self.env,
            self.service_labeling(),
        )

        self.log.info("Service started: %s", str(self.task_name))
        return self.cli.inspect_service(service["ID"])

    def cleanup(self):
        self.cli.remove_service(self.service["ID"])
        configs_json = self.service["Spec"]["TaskTemplate"]["ContainerSpec"].get(
            "Configs", []
        )
        for config_spec_json in configs_json:
            config_json = self.cli.inspect_config(config_spec_json["ConfigID"])
            config_labels = config_json["Spec"]["Labels"]

            if all(
                config_labels.get(key) == v
                for key, v in DockerSwarmOperator.shared_data_labeling().items()
            ):
                self.cli.remove_config(config_json["ID"])


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
def hacky_task_logs(cli: APIClient, task_id: str, since: int, is_tty: bool) -> bytes:
    params = {
        "details": False,
        "follow": False,
        "stdout": True,
        "stderr": True,
        "since": since,
        "timestamps": True,
        "tail": "all",
    }

    url = cli._url("/tasks/{0}/logs", task_id)
    res = cli._get(url, params=params, stream=False)
    res.raise_for_status()

    return cli._get_result_tty(False, res, is_tty)
