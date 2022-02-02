import base64
import io
import json
import os
import tarfile
from os.path import basename
from typing import Callable, Mapping
from typing import Optional, List, Any

from etna.operators.docker_operator import DockerOperator
from etna.operators.docker_operator_base import SwarmSharedData, DockerOperatorBase
from etna.operators.swarm_operator import DockerSwarmOperator

def run_on_docker(
        task_id: str,
        source_service: str,
        command: List[str],
        env: Mapping[str, str] = dict(),
        output_json: bool = False,
        output_b64: bool = False
) -> DockerOperatorBase:
    if os.environ.get('DEV_MODE'):
        return run_in_container(task_id, source_service, command, env=env, output_json=output_json, output_b64=output_b64)
    return run_in_swarm(task_id, source_service, command, env=env, output_json=output_json, output_b64=output_b64)


def run_in_swarm(
        task_id: str,
        source_service: str,
        command: List[str],
        env: Mapping[str, str] = dict(),
        include_external_network: bool = False,
        output_json: bool = False,
        output_b64: bool = False,
) -> DockerSwarmOperator:
    serialize_last_output = _prepare_input_outputs(
        output_b64, output_json
    )

    return DockerSwarmOperator(
        task_id=task_id,
        source_service=source_service,
        command=command,
        include_external_networks=include_external_network,
        serialize_last_output=serialize_last_output,
        env=env,
    )


def _prepare_input_outputs(output_b64, output_json):
    serialize_last_output: Optional[Callable[[bytes], Any]] = None
    if output_json:
        serialize_last_output = json.loads
    elif output_b64:
        serialize_last_output = base64.b64decode
    return serialize_last_output

def run_in_container(
        task_id: str,
        source_container: str,
        command: List[str],
        env: Mapping[str, str] = dict(),
        use_compose: bool = True,
        output_json: bool = False,
        output_b64: bool = False,
) -> DockerOperator:
    # Uses the containers created by the local development docker-compose system.
    if use_compose:
        source_container = f"monoetna_{source_container}_1"
    serialize_last_output = _prepare_input_outputs(
        output_b64, output_json
    )

    return DockerOperator(
        task_id=task_id,
        source_container_name=source_container,
        command=command,
        serialize_last_output=serialize_last_output,
        env=env,
    )
