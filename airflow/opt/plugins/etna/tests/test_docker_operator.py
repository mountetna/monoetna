import json
import logging
from dataclasses import dataclass
from io import StringIO
from unittest import mock

import pytest
from docker import APIClient

from etna import SwarmSharedData
from etna.operators.docker_operator import DockerOperator
from etna.tests.conftest import NotSoRandom


@dataclass
class FakeTaskInstance:
    task_id: str
    dag_id: str
    run_id: str


# In order to re-record these tests, you need to
# 1.  configure your local docker daemon into swarm mode (search this)
# 2.  run the test suite with a proxy from inside the container on port 8085 to the container's mounted docker socker
#        (actually, this is already setup.  use docker-compose-idea.yml with docker-compose, the mounted socket and
#         alternate entrypoint should work)
# 3.  Add the --record-mode argument to the pytest invocation (see pyvcr documentation)
@pytest.mark.vcr
@mock.patch("tempfile._Random", NotSoRandom)
def test_execute_docker_operator():
    cli = APIClient(base_url="http://localhost:8085")

    test_container_name = "test-container-operator"

    try:
        cli.remove_container(test_container_name)
    except:
        pass

    cli.create_container(image="bash", command=[], name=test_container_name)

    operator = DockerOperator(
        task_id="blahblahblah",
        source_container_name=test_container_name,
        command=[
            "bash",
            "-c",
            "[[ `cat /some/shared/data` == 'hello' ]] && for i in {0..10}; do echo $i; sleep 1; done",
        ],
        docker_base_url="http://localhost:8085",
        serialize_last_output=json.loads,
        swarm_shared_data=[
            SwarmSharedData(data=b"hello", remote_path="/some/shared/data")
        ],
    )

    test_buffer = StringIO()

    handler = logging.StreamHandler(test_buffer)
    operator.log.setLevel(logging.INFO)
    operator.log.addHandler(handler)

    context = dict(
        ti=FakeTaskInstance(
            task_id="task",
            dag_id="dag",
            run_id="run",
        )
    )
    operator.pre_execute(context)
    result = operator.execute(context)

    test_buffer.seek(0)
    assert test_buffer.read().split("\n") == [
        "Consuming logs now:",
        "0",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "Task terminated: exited",
        "Removing docker container: DockerOperator_a71d33f9732598817e8efb5693d985dd",
        "",
    ]

    assert result == 10
