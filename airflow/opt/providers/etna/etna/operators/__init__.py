import base64
import io
import json
import os
import tarfile
from os.path import basename
from typing import Callable, Mapping, Union, Type
from typing import Optional, List, Any

from serde.json import from_json

from .docker_operator import DockerOperator
from .docker_operator_base import DockerOperatorBase
from .swarm_operator import DockerSwarmOperator


def run_on_docker(
    task_id: str,
    source_service: str,
    command: List[str],
    env: Mapping[str, str] = dict(),
    output_json: Union[bool, Type] = False,
    output_b64: bool = False,
    docker_base_url="unix://var/run/docker.sock",
    include_external_networks=False,
    # NanoCpus (1 cpu = 1e9)
    cpu_limit: Optional[int] = int(1e9 * 0.5),
    # Bytes (1KB = 1024, 1MB = 1024 * 1024, 1GIG = 1024 * 1024 * 1024)
    mem_limit: Optional[int] = None,
    cpu_reservation: Optional[int] = None,
    mem_reservation: Optional[int] = None,
    allow_manager_nodes=False,
    additional_constraints: Optional[List[str]] = None,
) -> DockerOperatorBase:
    if os.environ.get("DEV_MODE"):
        return run_in_container(
            task_id,
            source_service,
            command,
            env=env,
            output_json=output_json,
            output_b64=output_b64,
            docker_base_url=docker_base_url,
        )
    return run_in_swarm(
        task_id,
        source_service,
        command,
        env=env,
        output_json=output_json,
        output_b64=output_b64,
        docker_base_url=docker_base_url,
        include_external_networks=include_external_networks,
        allow_manager_nodes=allow_manager_nodes,
        cpu_limit=cpu_limit,
        cpu_reservation=cpu_reservation,
        mem_reservation=mem_reservation,
        mem_limit=mem_limit,
        additional_constraints=additional_constraints
    )


def run_in_swarm(
    task_id: str,
    source_service: str,
    command: List[str],
    env: Mapping[str, str] = dict(),
    output_json: Union[bool, Type] = False,
    output_b64: bool = False,
    docker_base_url="unix://var/run/docker.sock",
    include_external_networks: Optional[bool] = False,
    # NanoCpus (1 cpu = 1e9)
    cpu_limit: Optional[int] = int(1e9 * 0.5),
    # Bytes (1KB = 1024, 1MB = 1024 * 1024, 1GIG = 1024 * 1024 * 1024)
    mem_limit: Optional[int] = None,
    cpu_reservation: Optional[int] = None,
    mem_reservation: Optional[int] = None,
    allow_manager_nodes=False,
    additional_constraints: Optional[List[str]] = None,
) -> DockerSwarmOperator:
    serialize_last_output = _prepare_input_outputs(output_b64, output_json)

    return DockerSwarmOperator(
        task_id=task_id,
        source_service=source_service,
        command=command,
        include_external_networks=include_external_networks,
        serialize_last_output=serialize_last_output,
        env=env,
        docker_base_url=docker_base_url,
        allow_manager_nodes=allow_manager_nodes,
        cpu_limit=cpu_limit,
        cpu_reservation=cpu_reservation,
        mem_reservation=mem_reservation,
        mem_limit=mem_limit,
        additional_constraints=additional_constraints
    )


def _prepare_input_outputs(output_b64, output_json):
    serialize_last_output: Optional[Callable[[bytes], Any]] = None
    if output_json:
        if isinstance(output_json, type):
            serialize_last_output = lambda bytes: from_json(
                output_json, bytes.encode("utf8")
            )
        else:
            serialize_last_output = json.loads
    elif output_b64:
        serialize_last_output = lambda bytes: bytes.decode("ascii")
    return serialize_last_output


def run_in_container(
    task_id: str,
    source_container: str,
    command: List[str],
    env: Mapping[str, str] = dict(),
    use_compose: bool = True,
    output_json: Union[bool, Type] = False,
    output_b64: bool = False,
    docker_base_url="unix://var/run/docker.sock",
) -> DockerOperator:
    # Uses the containers created by the local development docker-compose system.
    if use_compose:
        source_container = f"monoetna_{source_container}_1"
    serialize_last_output = _prepare_input_outputs(output_b64, output_json)

    return DockerOperator(
        task_id=task_id,
        source_container_name=source_container,
        command=command,
        serialize_last_output=serialize_last_output,
        env=env,
        docker_base_url=docker_base_url,
    )
