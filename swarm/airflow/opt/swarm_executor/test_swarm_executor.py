import datetime
import unittest
from typing import Optional, Callable

from airflow.utils.state import State

from . import SwarmExecutor
from airflow.models.taskinstance import TaskInstanceKey

class TestSwarmExecute(unittest.TestCase):
  def setUp(self) -> None:
    self.idx = 0
    self._teardown = []

  def tearDown(self) -> None:
    executor = self.create_executor()
    all_services = executor.executor_services()
    executor.log.info('Found %s services dangling', len(all_services))
    for service in all_services:
      executor.log.info('Cleaning up service %s', service.name)
      service.remove()

    for callable in self._teardown:
      callable()

  def on_tear_down(self, c: Callable):
    self._teardown.append(c)

  def new_id(self, template: str):
    self.idx += 1
    return template % self.idx

  def create_task_instance(self, dag_id: Optional[str] = None, task_id: Optional[str] = None, execution_date: Optional[datetime.datetime] = None, try_number = 1) -> TaskInstanceKey:
    return TaskInstanceKey(
      dag_id=dag_id or self.new_id('dag_%s'),
      task_id=task_id or self.new_id('task_%s'),
      execution_date=execution_date or datetime.datetime.now(),
      try_number=try_number,
    )

  def create_executor(self, job_id: Optional[str] = None):
    executor = SwarmExecutor()
    executor.job_id = job_id or self.new_id('executor_%s')
    executor.start()
    self.on_tear_down(lambda: executor.terminate())
    return executor

  def test_schedule_thing(self):
    executor1 = self.create_executor()
    executor2 = self.create_executor()
    task1 = self.create_task_instance()

    executor1.execute_async(
      task1,
      command=['airflow', 'tasks', 'run', 'true', 'some_parameter'],
      executor_config={'image': 'airflow:latest'}
    )
    executor1.sync()
    executor1.end()
    buffer = executor1.get_event_buffer()

    assert buffer == {task1: (State.SUCCESS, executor1.job_id)}