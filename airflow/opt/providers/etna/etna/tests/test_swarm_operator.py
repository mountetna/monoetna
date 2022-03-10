import json
import logging
from dataclasses import dataclass
from io import StringIO
from unittest import mock

import pytest
from docker import APIClient
from docker.types import (
    TaskTemplate,
    ContainerSpec,
    ConfigReference,
    SecretReference,
    Resources,
    Placement,
)

from etna.operators.docker_operator_base import (
    write_logs_and_yield_last,
    SwarmSharedData,
)
from etna.operators.swarm_operator import (
    create_service_definition,
    find_service,
    create_service_from_definition,
    find_local_network_ids,
    DockerSwarmOperator,
)


@dataclass
class FakeTaskInstance:
    task_id: str
    dag_id: str
    run_id: str


# @pytest.mark.vcr
# def test_swarm_cleanup():
#     test_service_name = "test-service-swarm"
#     test_config_name = "test-config"
#
#     cli = APIClient(base_url="http://localhost:8085")
#
#     try:
#         cli.remove_service(test_service_name)
#     except:
#         pass
#
#     try:
#         cli.remove_config(test_config_name)
#     except:
#         pass
#
#     config = cli.create_config(
#         test_config_name, b"abc", labels=DockerSwarmOperator.shared_data_labeling()
#     )
#     service = cli.create_service(
#         TaskTemplate(
#             ContainerSpec(
#                 image="alpine",
#                 command="true",
#                 configs=[ConfigReference(config["ID"], test_config_name)],
#             ),
#         ),
#         name=test_service_name,
#         labels=DockerSwarmOperator.service_labeling(),
#     )
#
#     while True:
#         tasks = cli.tasks(filters={"service": service["ID"]})
#         if not tasks:
#             continue
#         if tasks[0]["Status"]["State"] == "complete":
#             break
#
#     def get_configs():
#         return cli.configs(
#             filters={"label": join_labels(DockerSwarmOperator.shared_data_labeling())}
#         )
#
#     def get_services():
#         return cli.services(
#             filters={"label": join_labels(DockerSwarmOperator.service_labeling())}
#         )
#
#     assert len(get_configs()) > 0
#     assert len(get_services()) > 0
#
#     swarm_cleanup(cli)
#
#     assert len(get_configs()) > 0
#     assert len(get_services()) > 0
#
#     with freeze_time(datetime.datetime.now() + datetime.timedelta(hours=5)):
#         swarm_cleanup(cli)
#
#     assert len(get_configs()) == 0
#     assert len(get_services()) == 0


def test_write_logs_and_yield_last_with_intermittent_errors():
    def log_message(m: str, i: int) -> bytes:
        return f"2022-01-11T10:39:26.{str(i).zfill(9)}Z {m}".encode("utf8")

    def log_response(*args) -> bytes:
        return b"\n".join(args)

    response_iter = iter(
        [
            log_response(
                log_message("Hello!", 0),
                log_message("Yum yum", 1),
                log_message("Okay", 2),
            ),
            log_response(
                log_message("Okay", 1),
                b"Error grabbing logs: rpc error",
                log_message("Do not read me!", 5),
            ),
            log_response(
                log_message("Process done", 3),
            ),
        ]
    )

    def yield_response(*args):
        return next(response_iter)

    test_buffer = StringIO()
    handler = logging.StreamHandler(test_buffer)
    log = logging.getLogger("test-swarm-operator-logger")
    log.setLevel(logging.INFO)
    log.addHandler(handler)

    log_iter = write_logs_and_yield_last(yield_response, log)

    assert next(log_iter) == b"Okay"
    assert next(log_iter) == b"Process done"


# In order to re-record these tests, you need to
# 1.  configure your local docker daemon into swarm mode (search this)
# 2.  run the test suite with a proxy from inside the container on port 8085 to the container's mounted docker socker
#        (actually, this is already setup.  use docker-compose-idea.yml with docker-compose, the mounted socket and
#         alternate entrypoint should work)
# 3.  Add the --record-mode argument to the pytest invocation (see pyvcr documentation)
@pytest.mark.vcr
@mock.patch(
    "airflow.utils.strings.get_random_string", mock.MagicMock(return_value="abc")
)
def test_execute_swarm_operator():
    cli = APIClient(base_url="http://localhost:8085")
    # decorate_api_client(cli)

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
        task_id="blahblahblah",
        source_service=test_service_name,
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
    assert test_buffer.read().split("\n")[2:] == [
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
        "Task terminated: complete",
        "",
    ]

    assert result == 10


@pytest.mark.vcr
def test_find_and_create_service_definition():
    cli = APIClient(base_url="http://localhost:8085")

    test_service_name = "test-service-swarm"
    test_config_name = "test-config"
    test_secret_name = "test-secret"
    test_network_name = "test-network"
    test_network_name_global = "test-network-2"

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

    try:
        cli.remove_network(test_network_name_global)
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

    non_local_network = cli.create_network(
        test_network_name_global, driver="overlay", attachable=True
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
    create_service_from_definition(cli, definition, test_service_name, {}, {})
    cli.remove_service(test_service_name)

    assert definition.networks == [dict(Target=network["Id"])]

    cli.create_service(
        TaskTemplate(
            ContainerSpec(
                image="alpine",
            ),
            networks=[test_network_name_global],
        ),
        name=test_service_name,
    )
    service_data = find_service(cli, test_service_name)
    local_networks = find_local_network_ids(cli, service_data)
    assert local_networks == set()
    definition = create_service_definition(service_data, local_networks)

    assert definition.networks == []

    cli.remove_service(test_service_name)
    create_service_from_definition(cli, definition, test_service_name, {}, {})
    cli.remove_service(test_service_name)
