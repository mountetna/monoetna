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

from etna.operators.swarm_operator import (
    create_service_definition,
    find_service,
    create_service_from_definition,
    find_local_network_ids,
)


@pytest.mark.vcr
def test_find_and_create_service_definition():
    # cli = APIClient(base_url='http+unix://%2Fvar%2Frun%2Fdocker.sock/')
    # cli = APIClient(base_url='http+unix://%2Fvar%2Frun%2Fdocker.sock')
    # cli = APIClient(base_url="unix://var/run/docker.sock")
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
