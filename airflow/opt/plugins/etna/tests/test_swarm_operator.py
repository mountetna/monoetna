import logging
from dataclasses import dataclass
from io import StringIO

import json
import pytest
from airflow.models import TaskInstance
from docker import APIClient
from docker.types import (
    TaskTemplate,
    ContainerSpec,
    ConfigReference,
    SecretReference,
    Resources,
    Placement,
)

from etna.operators.swarm_operator import (
    create_service_definition,
    find_service,
    create_service_from_definition,
    find_local_network_ids, DockerSwarmOperator,
)


@dataclass
class FakeTaskInstance:
    task_id: str
    dag_id: str
    run_id: str

# In order to re-record these tests, you need to
# 1.  configure your lock docker daemon into swarm mode (search this)
# 2.  run the test suite with a proxy from inside the container on port 8085 to the container's mounted docker socker
#        (actually, this is already setup.  use docker-compose-idea.yml with docker-compose, the mounted socket and
#         alternate entrypoint should work)
# 3.  Add the --record-mode argument to the pytest invocation (see pyvcr documentation)
@pytest.mark.vcr
def test_execute_swarm_operator():
    cli = APIClient(base_url="http://localhost:8085")
    test_service_name = "test-service-swarm"

    try:
        cli.remove_service(test_service_name)
    except:
        pass

    # Try copying a service with many options set, and one with no options set, validating that we essentially 'copy'
    # their options e2e.
    cli.create_service(
        TaskTemplate(
            ContainerSpec(
                image="bash",
                command=[],
            ),
        ),
        name=test_service_name,
    )

    operator = DockerSwarmOperator(
        task_id='blahblahblah',
        source_service=test_service_name,
        command=["-c", "for i in {0..10}; do echo $i; sleep 1; done"],
        docker_base_url="http://localhost:8085",
        serialize_last_output=json.loads
    )

    test_buffer = StringIO()

    operator._log = logging.Logger('test', logging.INFO)
    operator.log.addHandler(logging.StreamHandler(test_buffer))

    operator.execute(dict(ti=FakeTaskInstance(
        task_id="task",
        dag_id="dag",
        run_id="run",
    )))

    assert test_buffer.read() == ""

@pytest.mark.vcr
def test_find_and_create_service_definition():
    cli = APIClient(base_url="http://localhost:8085")

    test_service_name = "test-service-swarm"
    test_config_name = "test-config"
    test_secret_name = "test-secret"
    test_network_name = "test-network"

    try:
        cli.remove_service(test_service_name)
    except:
        pass

    try:
        cli.remove_secret(test_secret_name)
    except:
        pass

    try:
        cli.remove_config(test_config_name)
    except:
        pass

    try:
        cli.remove_network(test_network_name)
    except:
        pass

    config = cli.create_config(test_config_name, b"abc")
    secret = cli.create_secret(test_secret_name, b"abc")
    network = cli.create_network(
        test_network_name,
        driver="overlay",
        attachable=True,
        labels={"com.docker.stack.namespace": "test"},
    )

    # Try copying a service with many options set, and one with no options set, validating that we essentially 'copy'
    # their options e2e.
    cli.create_service(
        TaskTemplate(
            ContainerSpec(
                image="alpine",
                command="true",
                user="0:0",
                mounts=["/usr/bin:/testmount:ro"],
                configs=[ConfigReference(config["ID"], test_config_name)],
                secrets=[SecretReference(secret["ID"], test_secret_name)],
                env=dict(test="abc"),
                tty=False,
            ),
            resources=Resources(mem_limit=5 * 1024 * 1024),
            networks=[test_network_name],
            placement=Placement(),
        ),
        name=test_service_name,
        labels={"com.docker.stack.namespace": "test"},
    )
    service_data = find_service(cli, test_service_name)
    definition = create_service_definition(
        service_data, find_local_network_ids(cli, service_data)
    )
    cli.remove_service(test_service_name)
    create_service_from_definition(cli, definition, test_service_name, {})
    cli.remove_service(test_service_name)

    assert definition.networks == [dict(Target=network["Id"])]

    cli.create_service(
        TaskTemplate(
            ContainerSpec(
                image="alpine",
            ),
        ),
        name=test_service_name,
    )
    service_data = find_service(cli, test_service_name)
    definition = create_service_definition(
        service_data, find_local_network_ids(cli, service_data)
    )

    assert definition.networks == []

    cli.remove_service(test_service_name)
    create_service_from_definition(cli, definition, test_service_name, {})
    cli.remove_service(test_service_name)
