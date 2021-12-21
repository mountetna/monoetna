import base64
import io
import json
import tarfile
from os.path import basename
from typing import Optional, List, Callable, Any

from etna.operators.docker_operator_base import SwarmSharedData
from etna.operators.swarm_operator import DockerSwarmOperator


def run_in_swarm(
    task_id: str,
    source_service: str,
    command: List[str],
    input_bytes: Optional[bytes] = None,
    input_file: Optional[str] = None,
    input_dir: Optional[str] = None,
    include_external_network: bool = False,
    output_json: bool = False,
    output_b64: bool = False,
):
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

    return DockerSwarmOperator(
        task_id=task_id,
        source_service=source_service,
        command=command,
        include_external_networks=include_external_network,
        swarm_shared_data=swarm_shared_data,
        serialize_last_output=serialize_last_output,
    )
