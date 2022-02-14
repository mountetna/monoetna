import functools
from datetime import datetime, timedelta
from typing import Union, Optional, Mapping, List

from airflow import DAG
from airflow.decorators import task
from airflow.operators.python import get_current_context

from etna.dags.decorators import dag
from etna.etls.etl_task_batching import enable_task_backfill, batch_end_context_key, get_batch_range, LOWEST_BOUND
from etna.etls.metis import load_metis_folders_batch, load_metis_files_batch
from etna.hooks.etna import EtnaHook, Folder, File
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

    interval = timedelta(hours=1)

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

            @task(do_xcom_push=False)
            def propagate_folder_updated_at(folders: List[Folder]):
                # For new etl startup, this step is not necessary.
                lower, _ = get_batch_range(get_current_context())
                if lower == LOWEST_BOUND:
                    return

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
                                   hook=hook, **inject_params))

        return etl(
            project_name=project_name,
            start_date=start_date,
            interval=interval,
            version=version,
        )(setup_tail)

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
