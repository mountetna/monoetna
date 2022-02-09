import functools
import os.path

from airflow.decorators import task
from pendulum import datetime
from datetime import timedelta
from typing import Mapping, Union, List
from typing import Optional

from airflow import DAG
from airflow.models.dag import (
    DagStateChangeCallback,
    ScheduleIntervalArgNotSet,
    ScheduleIntervalArg,
)

from etna.dags.callbacks import notify_slack_dag_callback
from etna.etls.batches import batch_end_context_key, enable_task_backfill
from etna.etls.metis import load_metis_folders_batch, load_metis_files_batch
from etna.hooks.etna import Folder, EtnaHook, File
from etna.utils.inject import inject
from etna.xcom.etna_xcom import pickled


def metis_etl(
        project_name: str,
        bucket_name: str,
        start_date: datetime,
        version: Union[int, str],
        propagate_updates=True,
        hook: Optional[EtnaHook] = None,
        inject_params: Mapping[str, str] = {},
):
    if hook is None:
        hook = EtnaHook.for_project(project_name)

    interval = timedelta(minutes=5)

    def instantiate_dag(fn):
        @functools.wraps(fn)
        def setup_tail():
            @task
            def tail_folders() -> List[Folder]:
                with hook.metis() as metis:
                    return pickled(load_metis_folders_batch(metis, bucket_name))

            @task
            def tail_files() -> List[File]:
                with hook.metis() as metis:
                    return pickled(load_metis_files_batch(metis, bucket_name))

            @task
            def propagate_folder_updated_at(folders: List[Folder]):
                with hook.metis() as metis:
                    for folder in folders:
                        children = metis.list_folder(project_name, bucket_name, folder.folder_path).folders
                        for child in children:
                            if child.updated_at_datetime < folder.updated_at_datetime:
                                metis.touch_folder(project_name, bucket_name, child.folder_path)

            folders = tail_folders()
            files = tail_files()

            if propagate_updates:
                propagate_folder_updated_at(folders)

            return inject(fn, dict(tail_folders=folders, tail_files=files, project_name=project_name, bucket_name=bucket_name,
                                   **inject_params))

        return etl(
            project_name=project_name,
            start_date=start_date,
            interval=interval,
            version=version,
        )(setup_tail)

    return instantiate_dag


def dag(
        on_failure_callback: Optional[DagStateChangeCallback] = None,
        on_success_callback: Optional[DagStateChangeCallback] = None,
        schedule_interval: ScheduleIntervalArg = ScheduleIntervalArgNotSet,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        inject_params: Mapping[str, str] = {},
        version: Union[int, str] = "",
        **kwds,
):
    """
    Creates a dag by running the wrapped function in a dag context whose dag_id will be the function's name,
    and whose description will be its doc string.
    :return: Dah'ya like dags?
    """

    def instantiate_dag(fn):
        with DAG(
                dag_id=fn.__name__ + str(version),
                description=fn.__doc__,
                on_failure_callback=on_failure_callback,
                on_success_callback=on_success_callback,
                start_date=start_date,
                end_date=end_date,
                schedule_interval=schedule_interval,
                **kwds,
        ) as dag:
            inject(fn, inject_params)
        return dag

    return instantiate_dag


def etl(project_name: str, start_date: datetime, interval: timedelta, version: Union[int, str],
        inject_params: Mapping = {}):

    def instantiate_dag(fn):
        nonlocal start_date

        # Stagger all etl start dates so as to reduce congestion.
        start_date += timedelta(seconds=abs(hash(fn.__name__)) % interval.seconds)

        new_dag: DAG = dag(
            start_date=start_date,
            schedule_interval=interval,
            catchup=True,
            inject_params=dict(
                batch_end_date="{{ " + batch_end_context_key + " }}",
                version=version,
                **inject_params
            ),
            default_args=dict(
                owner=project_name,
                retries=5,
            ),
            version=version,
        )(fn)

        for op in new_dag.tasks:
            enable_task_backfill(op)

        return new_dag

    return instantiate_dag


system_epoch = datetime(2021, 12, 22, 16, 56, 3, 185905)


# A dag context
def system_dag(interval: timedelta):
    def instantiate_dag(fn):
        return dag(
            start_date=system_epoch,
            schedule_interval=interval,
            default_args=dict(
                owner='administration',
                retries=5,
            ),
            # default_args=dict(
            on_failure_callback=notify_slack_dag_callback("failed: "),
            # ),
            catchup=False,
        )(fn)

    return instantiate_dag
