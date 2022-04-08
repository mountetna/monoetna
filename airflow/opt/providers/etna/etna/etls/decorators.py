import functools
import re
from datetime import timedelta
from typing import Union, Optional, Mapping, List

from airflow import DAG
from airflow.decorators import task
from airflow.operators.python import get_current_context

from etna.dags.decorators import dag, system_epoch
from etna.etls.etl_task_batching import (
    batch_end_context_key,
    get_batch_range,
    LOWEST_BOUND,
    batch_start_context_key,
)
from etna.etls.metis import (
    load_metis_folders_batch,
    load_metis_files_batch,
    MetisEtlHelpers,
)
from etna.hooks.etna import EtnaHook, Folder, File
from etna.utils.inject import inject
from etna.xcom.etna_xcom import pickled
from etna.etls.box import (
    BoxEtlHelpers,
    load_box_files_batch
)
from etna.hooks.box import BoxHook, Box


def metis_etl(
    project_name: str,
    bucket_name: str,
    version: Union[int, str],
    propagate_updates=True,
    hook: Optional[EtnaHook] = None,
    inject_params: Mapping[str, str] = {},
):
    """
    A decorator that converts a decorated function into a DAG by the same name, using all tasks instantiated within.

    This decorator is the main entry point for processing metis files with an etl.  The decorated function receives
    a number of useful XComArg objects and a helpers object to help set up metis tasks based on consuming changes
    scoped to the given project_name, bucket_name pair.

    1.  tail_folders -- an XComArg object that will resolve into a list of Folder objects inside tasks and link functions.
    2.  tail_files -- Similar to tail_folders, but with File objects.
    3.  helpers -- a MetisEtlHelpers object that contains useful functions for common etl processing.
    """
    if hook is None:
        hook = EtnaHook.for_project(project_name)

    interval = timedelta(days=1)

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
                result: List[str] = []
                # For new etl startup, this step is not necessary.
                lower, _ = get_batch_range(get_current_context())
                if lower == LOWEST_BOUND:
                    return []

                with hook.metis() as metis:
                    for folder in folders:
                        children = metis.list_folder(
                            project_name, bucket_name, folder.folder_path
                        ).folders
                        for child in children:
                            if child.updated_at_datetime < folder.updated_at_datetime:
                                result.append(child.folder_path)
                                metis.touch_folder(
                                    project_name, bucket_name, child.folder_path
                                )
                return result

            folders = tail_folders()
            files = tail_files()

            if propagate_updates:
                propagate_folder_updated_at(folders)

            helpers = MetisEtlHelpers(tail_folders=folders, tail_files=files, hook=hook)

            return inject(
                fn,
                dict(
                    tail_folders=folders,
                    tail_files=files,
                    project_name=project_name,
                    bucket_name=bucket_name,
                    helpers=helpers,
                    hook=hook,
                    **inject_params
                ),
            )

        return etl(
            project_name=project_name,
            interval=interval,
            version=version,
        )(setup_tail)

    return instantiate_dag


def box_etl(
    folder_name: str,
    version: Union[int, str],
    hook: Optional[BoxHook] = None,
    project_name: str = "administration",
    inject_params: Mapping[str, str] = {},
):
    """
    A decorator that converts a decorated function into a DAG by the same name, using all tasks instantiated within.

    This decorator is the main entry point for processing Box files with an etl.  The decorated function receives
    a number of useful XComArg objects and a helpers object to help set up Box tasks based on consuming changes
    scoped to the given folder in Box.

    1.  tail_files -- an XComArg object that will resolve into a list of File objects inside tasks.
    """
    if not Box.valid_folder_name(folder_name):
        raise ValueError(f"Invalid folder name: {folder_name}. Only alphanumeric characters, -, and spaces are allowed.")

    if hook is None:
        hook = BoxHook.for_project(project_name)

    interval = timedelta(days=1)

    def instantiate_dag(fn):
        @functools.wraps(fn)
        def setup_tail():
            @task
            def tail_files() -> List[File]:
                with hook.box() as box:
                    return pickled(load_box_files_batch(box, folder_name))

            files = tail_files()

            helpers = BoxEtlHelpers(
                hook=hook
            )

            return inject(
                fn,
                dict(
                    tail_files=files,
                    helpers=helpers,
                    **inject_params
                ),
            )

        return etl(
            project_name=project_name,
            interval=interval,
            version=version,
        )(setup_tail)

    return instantiate_dag


def etl(
    project_name: str,
    interval: timedelta,
    version: Union[int, str],
    inject_params: Mapping = {},
):
    def instantiate_dag(fn):
        start_date = system_epoch
        # Stagger start times for etls
        start_date += timedelta(
            seconds=hash(fn.__name__) % int(interval.total_seconds())
        )

        new_dag: DAG = dag(
            start_date=start_date,
            schedule_interval=interval,
            catchup=False,
            inject_params=dict(
                batch_start_date="{{ " + batch_start_context_key + " }}",
                batch_end_date="{{ " + batch_end_context_key + " }}",
                version=version,
                **inject_params
            ),
            default_args=dict(
                owner=project_name,
                depends_on_past=True,
                retries=3,
            ),
            tags=[project_name, 'etl'],
            version=version,
        )(fn)

        return new_dag

    return instantiate_dag
