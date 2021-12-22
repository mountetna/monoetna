import tempfile
from typing import Optional, Callable, Tuple, Set, Mapping

import docker.errors
from airflow import AirflowException
from docker import APIClient
from docker.models.containers import Container

from etna.operators.docker_operator_base import DockerOperatorBase


def find_container(cli: APIClient, container_name: str) -> Optional[Container]:
    try:
        return Container(cli.inspect_container(container_name))
    except docker.errors.APIError as e:
        if e.response.status_code != 404:
            raise e
    return None


class DockerOperator(DockerOperatorBase):
    source_container_name: str
    source: Container

    def __init__(
        self,
        *args,
        source_container_name: str,
        **kwds,
    ):
        self.source_container_name = source_container_name
        super(DockerOperator, self).__init__(*args, **kwds)

    def pre_execute(self, context) -> None:
        super(DockerOperator, self).pre_execute(context)
        self.source = find_container(self.cli, self.source_container_name)
        if not self.source:
            raise AirflowException(f"Could not find source container: {self.source_container_name}")

    @staticmethod
    def container_labeling():
        return {"etna.operators.docker_operator.container": "true"}

    def cleanup(self) -> None:
        if self.cli is not None:
            self.log.info("Removing docker container: %s", self.task_name)
            self.cli.remove_container(self.task_name)

    def _find_or_create_container(self) -> Mapping:
        existing = find_container(self.cli, self.task_name)
        if existing is not None:
            return existing.attrs

        config = self.source.attrs.get("Config", {})

        new_host_config = {
            k: v
            for k, v in self.source.attrs["HostConfig"].items()
            if k not in {"PortBindings", "RestartPolicy", "AutoRemove", "LogConfig"}
        }
        binds = new_host_config.setdefault("Binds", [])
        if binds is None:
            binds = []
            new_host_config['Binds'] = binds

        for data in self.swarm_shared_data:
            with tempfile.NamedTemporaryFile(delete=False) as file:
                file.write(data.data)

            binds.append(f"{file.name}:{data.remote_path}:ro")

        container = self.cli.create_container(
            self.source.attrs.get('ImageID', self.source.attrs['Image']),
            self.command,
            entrypoint=config.get('Entrypoint', []),
            user=config.get("User") or None,
            tty=config.get("Tty"),
            environment=config.get("Env", []) + list(f"{k}={v}" for k, v in self.env.items()),
            host_config=new_host_config,
            labels=self.container_labeling(),
            name=self.task_name,
        )
        for network_name, network in (
                self.source.attrs["NetworkSettings"].get("Networks", {}).items()
        ):
            self.cli.connect_container_to_network(
                container["Id"],
                network["NetworkID"] or network_name,
                links=network.get("Links", None),
                driver_opt=network.get("DriverOpt", None),
            )

        return container

    def _start_task(self):
        container = self._find_or_create_container()
        self.log.info("Starting %s inside of %s", self.command, self.source_container_name)
        self.cli.start(container["Id"])

    def _check_task(self) -> Tuple[str, Optional[str]]:
        container = self.cli.inspect_container(self.task_name)
        state = container["State"]
        if state.get("ExitCode"):
            return "failed", state.get("Error")

        return state["Status"], state.get("Error")

    @property
    def _next_log_batch_since(self) -> Callable[[int], bytes]:
        def next_log_batch(since: int):
            return hacky_container_logs(self.cli, self.task_name, since)

        return next_log_batch

    @property
    def successful_states(self) -> Set[str]:
        return {
            "exited",
        }

    @property
    def completed_states(self) -> Set[str]:
        return self.successful_states.union({"failed", "dead", "paused"})


# Addresses issue with current python client's container logs method.
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
def hacky_container_logs(cli: APIClient, container: str, since: int) -> bytes:
    params = {
        "stderr": 1,
        "stdout": 1,
        "timestamps": 1,
        "follow": 0,
        "tail": "all",
    }

    if since > 0:
        params["since"] = since
    else:
        params["since"] = 1

    url = cli._url("/containers/{0}/logs", container)
    res = cli._get(url, params=params, stream=False)
    return cli._get_result(container, False, res)
