import time
from dataclasses import dataclass
from typing import Optional, List, Dict

import docker.errors
from airflow.models import TaskInstance
from airflow.utils.strings import get_random_string
from docker import types, APIClient

from airflow.providers.docker.operators.docker_swarm import DockerSwarmOperator as OrigDockerSwarmOperator
from docker.models.services import Service
from docker.types import ConfigReference

class SwarmSharedData(dataclass):
  data: bytes
  remote_path: str

def join_labels(labels: Dict[str, str]) -> List[str]:
  return [ f"{k}=={v}" for k, v in labels.items() ]

# TODO: Executor can hook this up.
def swarm_cleanup(cli: APIClient):
  for config in cli.configs(filters={'label': join_labels(DockerSwarmOperator.shared_data_labeling())}):
    try:
      config.remove()
    except  docker.errors.APIError:
      pass  # ignore, maybe being used still.

  for service in cli.services(filters={'label': join_labels(DockerSwarmOperator.service_labeling())}):
    task_status = cli.tasks(filters={'service': service['ID']})[0]['Status']
    if task_status['State'] in ['complete', 'failed', 'shutdown', 'rejected', 'orphaned', 'remove']:
      # timestamp: str = task_status['Timestamp']
      # timestamp = datetime.strptime(timestamp, datetime.isoformat())
      # datetime.now() - datetime.
      # TODO: Maybe remove based on timestamp?  Also retry around APIErrors, especially version errors
      service.remove()

    pass

class DockerSwarmOperator(OrigDockerSwarmOperator):
  swarm_shared_data: List[SwarmSharedData]
  cli: APIClient
  ti: TaskInstance

  def __init__(self, *args, swarm_shared_data: Optional[List[SwarmSharedData]], **kwds):
    super(DockerSwarmOperator, self).__init__(*args, **kwds)
    self.swarm_shared_data = swarm_shared_data or []

  def execute(self, context) -> None:
    self.cli = self._get_cli()
    self.environment['AIRFLOW_TMP_DIR'] = self.tmp_dir
    self.ti = context['ti']
    return self._run_service()

  @staticmethod
  def shared_data_labeling():
    return {'etna.operators.swarm_operator.shared_data': 'true'}

  @staticmethod
  def service_labeling():
    return {'etna.operators.swarm_operator.service': 'true'}

  def _service_name(self):
    return f"{self.ti.dag_id}-{self.ti.task_id}-{self.ti.run_id}"

  def _prepare_shared_data(self):
    configs: List[ConfigReference] = []

    for shared_data in self.swarm_shared_data:
      config = self.cli.create_config(
        f"{self.dag_id}-{self.task_id}-shared-data-{get_random_string()}",
        shared_data.data,
        self.shared_data_labeling(),
      )

      configs.append(ConfigReference(
        config.id,
        config.name,
        filename=shared_data.remote_path
      ))

    return configs

  def _find_or_create_service(self) -> Service:
    for service in self.cli.services(filters={'name': self._service_name()}):
      return service

    service = self.cli.create_service(
      types.TaskTemplate(
        container_spec=types.ContainerSpec(
          image=self.image,
          command=self.format_command(self.command),
          mounts=self.mounts,
          env=self.environment,
          user=self.user,
          tty=self.tty,
          configs=self.configs + self._prepare_shared_data(),
          secrets=self.secrets,
        ),
        restart_policy=types.RestartPolicy(condition='none'),
        resources=types.Resources(mem_limit=self.mem_limit),
        networks=self.networks,
        placement=self.placement,
      ),
      name=self._service_name(),
      labels=self.service_labeling(),
      mode=self.mode,
    )

    self.log.info('Service started: %s', str(self._service_name()))

    return service

  def on_kill(self) -> None:
    if self.cli is not None:
      self.log.info('Removing docker service: %s', self.service['ID'])
      self.cli.remove_service(self.service['ID'])

  def _run_service(self) -> None:
    self.log.info('Starting docker service from image %s', self.image)
    if not self.cli:
      raise Exception("The 'cli' should be initialized before!")

    self.service = self._find_or_create_service()

    # wait for the service to start the task
    while not self.cli.tasks(filters={'service': self.service['ID']}):
      time.sleep(1)

    if self.enable_logging:
      self._stream_logs_to_output()

    while True:
      if self._has_service_terminated():
        self.log.info('Service status before exiting: %s', self._service_status())
        break

    if self.service and self._service_status() != 'complete':
      if self.auto_remove:
        self.cli.remove_service(self.service['ID'])
      raise AirflowException('Service did not complete: ' + repr(self.service))
    elif self.auto_remove:
      if not self.service:
        raise Exception("The 'service' should be initialized before!")
      self.cli.remove_service(self.service['ID'])