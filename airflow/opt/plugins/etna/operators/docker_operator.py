import tempfile
from typing import Optional, Callable, Any, Tuple, Set

from docker import APIClient
from docker.models.containers import Container

from etna.operators.docker_operator_base import DockerOperatorBase


def find_container(cli: APIClient, container_name: str) -> Optional[Container]:
  containers = cli.containers(limit=1, filters=dict(name=container_name))
  if containers:
    return Container(attrs=containers[0])
  return None

class DockerSwarmOperator(DockerOperatorBase):
  source_container_name: str
  tty: bool
  source: Container

  def __init__(
          self,
          *args,
          source_container_name: str,
          **kwds,
  ):
    self.source_container_name = source_container_name
    super(DockerSwarmOperator, self).__init__(*args, **kwds)

  def pre_execute(self, context) -> None:
    super(DockerSwarmOperator, self).pre_execute(context)
    self.source = find_container(self.cli, self.source_container_name)
    self.tty = self.source.attrs['Config']['Tty']

  @staticmethod
  def container_labeling():
    return {"etna.operators.docker_operator.container": "true"}

  def cleanup(self) -> None:
    if self.cli is not None:
      self.log.info("Removing docker container: %s", self.task_name)
      self.cli.remove_container(self.task_name)


  def execute(self, context: Any):
    config = self.source.attrs['Config']

    new_host_config = {k: v for k, v in self.source.attrs['HostConfig'] if k not in {'PortBindings', 'RestartPolicy', 'AutoRemove', 'LogConfig'}}
    binds = new_host_config.setdefault('Binds', [])

    for data in self.swarm_shared_data:
      with tempfile.NamedTemporaryFile(delete=False) as file:
        file.write(data.data)

      binds.append(f"{file.name}:{data.remote_path}:ro")

    container = self.cli.create_container(
      self.source.image,
      self.command,
      user=config['User'] or None,
      tty=self.tty,
      environment=config['Env'] + list(f"{k}={v}" for k, v in self.env.items()),
      host_config=new_host_config,
      labels=self.container_labeling(),
      name=self.task_name,
    )

    for network_name, network in self.source.attrs['NetworkSettings'].get('Networks', {}).items():
      self.cli.connect_container_to_network(
        container['Id'],
        network['NetworkID'],
        links=network.get('Links', None),
        driver_opt=network.get('DriverOpt', None),
      )

    return self._run()

  def _check_task(self) -> Tuple[str, Optional[str]]:
    containers = self.cli.containers(limit=1, filters=dict(name=self.task_name))
    if containers:
      state = containers[0]['State']
      if state.get('ExitCode'):
        return 'failed', state.get('Error')

      return state['Status'], state.get('Error')

  def _next_log_batch_since(self) -> Callable[[int], bytes]:
    def next_log_batch(since: int):
      return hacky_container_logs(self.cli, self.task_name, since)

    return next_log_batch

  @property
  def successful_states(self) -> Set[str]:
    return {"exited",}

  @property
  def completed_states(self) -> Set[str]:
    return self.successful_states.union({"failed", "dead", "paused"})

# Addresses issue with current python client's container logs method.
# We force the entire logs to be consumed in a singular buffer rather than streamed.
# This isn't ideal but the streaming solution offered in the underlying implementation
# has many problems.  This should be ok since we still grab chunks using the 'since'
# property and are capable of consuming logs "as they are produced".  Worst case scenario,
# an extremely prolific logging service might cause network and memory bloating with
# the buffered approach, but this would be indicative of a deeper problem.
# In future maybe we can rewrite the python docker library's streaming implementation
# to fix underlying issues:
#   1. ability to cancel streaming (currently it removes socket timeout and does not allow retrieval
#      of underlying socket to cancel it, resulting in a leaking stream when services are completed)
#   2. it attempts to hack into internals of requests implementation to retrieve the underlying socket,
#      which is not portable when using the vcr requests library, nor is it portable when using an http
#      adapters or even different python versions.  This is a flaw in the way python http is implemented,
#      its attempt to hide underlying sockets causes grief (ie: stop hiding resources that are not actually private)
def hacky_container_logs(cli: APIClient, container: str, since: int) -> bytes:
  params = {
    'stderr': 1,
    'stdout': 1,
    'timestamps': 1,
    'follow': 0,
    'tail': 'all',
  }

  if since > 0:
    params['since'] = since

  url = cli._url("/containers/{0}/logs", container)
  res = cli._get(url, params=params, stream=False)
  return cli._get_result(container, False, res)
