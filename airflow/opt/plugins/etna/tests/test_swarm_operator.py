import pytest
from docker import APIClient
from docker.types import TaskTemplate, ContainerSpec

from etna.operators.swarm_operator import create_service_definition, find_service


@pytest.mark.vcr
def test_find_and_create_service_definition():
  cli = APIClient(base_url='http+unix://%2Fvar%2Frun%2Fdocker.sock/')

  test_service_name = "test-service-swarm"

  try:
    cli.remove_service(test_service_name)
  except:
    pass

  service = cli.create_service(TaskTemplate(
    ContainerSpec(
      image='alpine',
      command='true',
    ),
  ), name=test_service_name)

  try:
    service_definition = create_service_definition(find_service(cli, test_service_name))
  finally:
    service.remove()
