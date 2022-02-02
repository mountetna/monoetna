import base64
import io
import json
import tarfile
from os.path import basename
from typing import Callable, Mapping
from typing import Optional, List, Any

from etna.operators.docker_operator import DockerOperator
from etna.operators.docker_operator_base import SwarmSharedData
from etna.operators.swarm_operator import DockerSwarmOperator


def run_in_swarm(
        task_id: str,
        source_service: str,
        command: List[str],
        env: Mapping[str, str] = dict(),
        input_bytes: Optional[bytes] = None,
        input_file: Optional[str] = None,
        input_dir: Optional[str] = None,
        include_external_network: bool = False,
        output_json: bool = False,
        output_b64: bool = False,
) -> DockerSwarmOperator:
    serialize_last_output, swarm_shared_data = _prepare_input_outputs(
        input_bytes, input_dir, input_file, output_b64, output_json
    )

    return DockerSwarmOperator(
        task_id=task_id,
        source_service=source_service,
        command=command,
        include_external_networks=include_external_network,
        swarm_shared_data=swarm_shared_data,
        serialize_last_output=serialize_last_output,
        env=env,
    )


def _prepare_input_outputs(input_bytes, input_dir, input_file, output_b64, output_json):
    swarm_shared_data: List[SwarmSharedData] = []
    serialize_last_output: Optional[Callable[[bytes], Any]] = None
    if input_bytes:
        swarm_shared_data.append(
            SwarmSharedData(
                data=input_bytes,
                remote_path="/input",
            )
        )
    elif input_file:
        swarm_shared_data.append(
            SwarmSharedData(
                data=open(input_file, "r").read(),
                remote_path="/input/" + basename(input_file),
            )
        )
    elif input_dir:
        buff = io.BytesIO(b"")
        tar = tarfile.open(fileobj=buff)
        tar.add(input_dir, arcname=basename(input_dir))
        tar.close()

        swarm_shared_data.append(
            SwarmSharedData(
                data=buff.read(),
                remote_path="/input/" + basename(input_dir) + ".tar",
            )
        )
    if output_json:
        serialize_last_output = json.loads
    elif output_b64:
        serialize_last_output = base64.b64decode
    return serialize_last_output, swarm_shared_data

def run_in_container(
        task_id: str,
        source_container: str,
        command: List[str],
        env: Mapping[str, str] = dict(),
        input_bytes: Optional[bytes] = None,
        use_compose: bool = True,
        input_file: Optional[str] = None,
        input_dir: Optional[str] = None,
        output_json: bool = False,
        output_b64: bool = False,
) -> DockerOperator:
    # Uses the containers created by the local development docker-compose system.
    if use_compose:
        source_container = f"monoetna_{source_container}_1"
    serialize_last_output, swarm_shared_data = _prepare_input_outputs(
        input_bytes, input_dir, input_file, output_b64, output_json
    )

    return DockerOperator(
        task_id=task_id,
        source_container_name=source_container,
        command=command,
        swarm_shared_data=swarm_shared_data,
        serialize_last_output=serialize_last_output,
        env=env,
    )
