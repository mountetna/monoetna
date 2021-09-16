# https://airflow.apache.org/docs/apache-airflow/stable/executor/index.html
# https://github.com/apache/airflow/blob/main/airflow/executors/kubernetes_executor.py
import base64
import functools
import hashlib
import itertools
import json
import os
import pickle
import time
from typing import Any, List, Optional, NamedTuple, Mapping, Callable, Generator, \
  MutableMapping, Iterator, Tuple
from airflow.configuration import conf

import docker
import docker.errors
import docker.models
import docker.types
from airflow.exceptions import AirflowException
from airflow.executors.base_executor import BaseExecutor, CommandType
from airflow.models.taskinstance import TaskInstance, TaskInstanceKey
from airflow.utils.event_scheduler import EventScheduler
from airflow.utils.log.logging_mixin import LoggingMixin
from airflow.utils.state import State
from docker.models.services import Service
from docker.types import RestartPolicy
from docker.types.services import RestartConditionTypesEnum


class QueuedSwarmCommand(NamedTuple):
  key: TaskInstanceKey
  command: CommandType
  queue: Optional[str]
  executor_config: Optional[Mapping[str, Any]]

  @property
  def as_service_name(self):
    # A unique key that will ensure only one creation of a task instance
    # in the swarm, in any case.
    parts = [
      "dag",
      self.key.dag_id,
      "task",
      self.key.task_id,
      str(self.key.execution_date.timestamp()),
      str(self.key.try_number)
    ]

    # Ensures safe encoding
    return hashlib.md5(':'.join(parts).encode('utf-8')).hexdigest()


def task_instance_key_from_service(service: Service) -> TaskInstanceKey:
  labels = service.attrs['Spec']['Labels']
  task_id = labels.get('airflow_executor.task_id', '')
  dag_id = labels.get('airflow_executor.dag_id', '')
  execution_date = labels.get('airflow_executor.execution_date', '')
  if execution_date:
    execution_date = pickle.loads(base64.b64decode(execution_date))

  try_number = int(labels.get('airflow_executor.try_number', 1))

  return TaskInstanceKey(
    dag_id=dag_id,
    task_id=task_id,
    execution_date=execution_date,
    try_number=try_number,
  )


def task_instance_key_as_labels(key: TaskInstanceKey) -> Mapping[str, str]:
  return {
    "airflow_executor.task_id": key.task_id,
    "airflow_executor.dag_id": key.dag_id,
    "airflow_executor.execution_date":
      base64.b64encode(pickle.dumps(key.execution_date)).decode('utf-8'),
    "airflow_executor.try_number": str(key.try_number),
  }


def service_matches_label(service: Service, label: Mapping[str, str]):
  spec_labels = service.attrs['Spec']['Labels']
  return all(spec_labels.get(k) == v for k, v in label.items())


failed_stats = {'failed', 'rejected', 'orphaned'}


def service_status(service: Service) -> Optional[bool]:
  """
  :param service:
  :return: True if the service completed successfully, False if it failed,
  None if it is running.
  """
  all_complete = None
  for task in service.tasks():
    task_state = task['Status']['State']
    if task_state == 'complete':
      if all_complete is None:
        all_complete = True
      continue

    all_complete = False
    if task_state in failed_stats:
      return False

  if all_complete:
    return True

  return None


def label_as_string(label: Mapping[str, str]):
  for k, v in label.items():
    return f"{k}={v}"
  return ""


class SwarmExecutor(BaseExecutor, LoggingMixin):
  """
  Executor for DockerSwarm.
  Each airflow task is mapped to a unique, deterministic service in docker swarm.
    Within that service is a single docker task whose status is synced with
    the airflow scheduler's status.
  """

  def __init__(self, heartbeat_rate: Optional[int] = None):
    self.event_scheduler: Optional[EventScheduler] = None
    self.docker_client: docker.DockerClient = docker.DockerClient(
      base_url='unix://var/run/docker.sock'
    )
    self.heartbeat_rate = heartbeat_rate or conf.getint('scheduler', 'SCHEDULER_HEARTBEAT_SEC')
    self.new_service_queue: List[QueuedSwarmCommand] = []
    self.queue_configs: MutableMapping[str, Any] = {}
    super().__init__()

  def get_queue_config(self, queue: str):
    if queue not in self.queue_configs:
      path = os.path.realpath(f"/run/swarm_executor/{queue}.json")
      if not path.startswith('/run/swarm_executor'):
        raise AirflowException('Invalid queue name: ' + queue)
      if os.path.isfile(path):
        self.queue_configs[queue] = json.load(open(path, 'r'))

    return self.queue_configs.get(queue, {})

  def has_task(self, task_instance: TaskInstance) -> bool:
    return super().has_task(task_instance) or task_instance.queued_by_job_id == self.job_id

  def _reap_finished_services(self, services: List[Service]) -> Iterator[Callable]:
    for service in services:
      status = service_status(service)
      if status is None:
        continue

      target_state: str = State.SUCCESS if status is True else State.FAILED
      yield functools.partial(self._reap_finished_service, target_state, service)

  def _reap_finished_service(self, target_state: str, service: Service):
    key = task_instance_key_from_service(service)
    stderr = b'\n'.join(service.logs(stderr=True, stdout=True))
    if self.remove_service(service):
      self.log.info('Swarm service completed, reaping %s', key)
      if stderr.strip():
        try:
          readable = stderr.decode('utf-8')
          self.log.error('Stderr from completed job:\n%s', readable)
        except UnicodeDecodeError:
          self.log.error('Stderr from completed job:\n%s', stderr)
        self.change_state(key, target_state, self.job_id)
    # If by chance we actually have this key, we've now discovered it has completed,
    # so our internal service queue should clear it.
    # This might happen if a create service command fails (connection error) but
    # in reality the service starts and is completed.
    queued = [q for q in self.new_service_queue if q.key == key]
    for q in queued:
      self.new_service_queue.remove(q)

  def _start_new_services(self, services: List[Service]):
    for new_service_item in self.new_service_queue[:]:
      config: MutableMapping[str, Any] = {}
      if new_service_item.executor_config:
        config.update(new_service_item.executor_config)
      config.update(self.get_queue_config('default'))
      config.update(self.get_queue_config(new_service_item.queue))

      image = config.pop('image')
      command = new_service_item.command

      config.setdefault('labels', {})
      config['labels'].update(
        **self.executor_label(),
        **self.scheduler_label(self.job_id),
        **task_instance_key_as_labels(new_service_item.key),
      )

      config['name'] = new_service_item.as_service_name

      yield functools.partial(self._start_service,
                              new_service_item, command, image, config)

  def _start_service(self,
                     new_service_item: QueuedSwarmCommand,
                     command: CommandType,
                     image: str,
                     config: Mapping[str, Any]):
    try:
      if command[0:3] != ["airflow", "tasks", "run"]:
        raise ValueError('The command must start with ["airflow", "tasks", "run"].')

      self.docker_client.configs.create()

      self.docker_client.services.create(
        image,
        command,
        maxreplicas=1,
        restart_policy=RestartPolicy(
          condition=RestartConditionTypesEnum.NONE,
        ),
        **config
      )
    except docker.errors.APIError as e:
      # this could happen if the service was started, but our connection
      # fails, for instance.  In that case, it may still run and complete,
      # so we wouldn't want to mark it as failure but then switch to success
      # later on.  Instead, we let the service exist and let it converge
      # naturally with a future reap.
      if e.response.status_code != 409:
        self.log.error('Failed to start up new service', exc_info=True)
        # we should only mark failure when we can explicitly confirm that
        # creating the service has failed.
        self.fail(new_service_item.key, e)
    except (ValueError, TypeError) as e:
      self.log.error('Failed with service configuration', exc_info=True)
      self.fail(new_service_item.key, e)
    else:
      self.log.info('New service is running: %s', new_service_item.key)
      self.running.add(new_service_item.key)

    self.log.info('Processed starting new task: %s', new_service_item.key)
    self.new_service_queue.remove(new_service_item)

  def _wait_for_pending_services(self, services: List[Service], await_pending=False):
    def wait():
      time.sleep(1)

    if services and await_pending:
      yield wait

  def work_to_converge(self, await_pending = False):
    while True:
      services = self.executor_services(filter_label=self.scheduler_label(self.job_id))
      self.log.info('Found %s services belonging to this scheduler', len(services))

      found_some = False
      for unit in itertools.chain(
              self._reap_finished_services(services),
              self._start_new_services(services),
              self._wait_for_pending_services(services, await_pending),
      ):
        found_some = True
        yield unit

      if not found_some:
        break

  def start(self) -> None:
    if not self.job_id:
      raise AirflowException("Could not get scheduler_job_id")
    self.event_scheduler = EventScheduler()

  def execute_async(
          self,
          key: TaskInstanceKey,
          command: CommandType,
          queue: Optional[str] = None,
          executor_config: Optional[Any] = None,
  ) -> None:
    """Executes task asynchronously"""
    self.log.info('Add task %s with command %s with executor_config %s', key, command, executor_config)
    self.event_buffer[key] = (State.QUEUED, self.job_id)
    self.new_service_queue.append(QueuedSwarmCommand(
      key=key,
      command=command,
      queue=queue,
      executor_config=executor_config,
    ))

  def sync(self) -> None:
    """Synchronize task state."""
    start = time.time()
    for unit in self.work_to_converge():
      # Keep processing and waiting for pending work up until we're approaching
      # our heartbeat rate.
      if time.time() - start > self.heartbeat_rate:
        break

      try:
        unit()
      except:
        self.log.error('Failed to process some work', exc_info=True)
    self.event_scheduler.run(blocking=False)

  def executor_services(self,
                        filter_label: Optional[Mapping[str, str]] = None) -> List[Service]:
    executor_services = self.docker_client.services.list(
      filters=dict(label=label_as_string(self.executor_label())))

    if filter_label:
      return [s for s in executor_services
              if service_matches_label(s, filter_label)]

    return executor_services

  def try_adopt_task_instances(self, tis: List[TaskInstance]) -> List[TaskInstance]:
    """
    Called by airflow's SchedulerJob periodically with tasks that do not have a
    owning scheduler job id, either because it is none or the state of that job
    is now not running.  The task instances themselves are either queued or running.
    """
    adoptable_task_instances = {ti.key: ti for ti in tis if ti.queued_by_job_id}
    adopted_task_instances = {ti.key for ti in adoptable_task_instances}

    to_reset = []  # these will be task instances that we do not adopt, and should be reset.
    # First identify those task instances we are not going to adopt, warn about them,
    # and add them to the reset group.
    for task_instance in tis:
      # These are apparently historical tasks that did not record their job id, so we
      # should ignore these in general and let the resetting work handle them.
      if task_instance.key not in adoptable_task_instances:
        self.log.warning('Found task_instance with null queued_by_job_id in try_adopt_task_instances: %s',
                         task_instance)
        to_reset.append(task_instance)

    # Keep converging all services until they are adopted into loving scheduler jobs.
    # Basically, it's possible for the airflow db's sense of scheduler job ownership
    # to not match the actual state of the swarm.  So we search all executor services
    # for services that could be adopted according to the airflow db and formalize that
    # adoption, even if it has already occurred partially (swarm only) elsewhere.
    services = self.executor_services()

    for service in services:
      ti_key = task_instance_key_from_service(service)
      if ti_key in adoptable_task_instances:
        if self.adopt_service(service):
          adopted_task_instances.add(ti_key)
          self.running.add(ti_key)

    # Return any of the original task instances we did not, actually, adopt.
    return [
      ti for ti in tis
      if ti.key not in adopted_task_instances
    ]

  def remove_service(self, service: Service) -> bool:
    try:
      if not service_matches_label(service, self.scheduler_label(self.job_id)):
        return False
      service.remove()
    except docker.errors.APIError as e:
      if e.response.status_code in {409, 404}:
        return False
      raise
    return True

  def adopt_service(self, service: Service) -> bool:
    """
    Attempts to 'adopt' a docker swarm service to this scheduler, that is
    by labeling it with the unique job id of this scheduler.
    This can fail with concurrent updates due to the optimistic locking mechanism.
    In that case, False is returned.
    Otherwise, True is returned on successful adoption.
    :param service:
    :return:
    """
    try:
      service.update(labels=self.scheduler_label(self.job_id))
    except docker.errors.APIError as e:
      if e.response.status_code == 409:
        service.reload()
        return False
      raise

    return True

  def scheduler_label(self, job_id=None):
    return {'airflow_executor.scheduler_job': job_id or self.job_id}

  def executor_label(self):
    return {'airflow_executor': 'SwarmExecutor1'}

  def end(self) -> None:
    """Synchronously complete all work"""
    try:
      for unit in self.work_to_converge(True):
        unit()
    finally:
      self.docker_client.close()

  def terminate(self):
    """Terminate the executor is not doing anything."""
    self.docker_client.close()
