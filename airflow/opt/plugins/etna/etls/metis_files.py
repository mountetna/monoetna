import functools
from datetime import timedelta, datetime, timezone
from typing import Union, Literal, List, Callable, Tuple

from airflow import DAG
from airflow.decorators import task
from airflow.models.taskinstance import Context
from airflow.operators.python import get_current_context

from etna.dags import etl, run_in_swarm
from etna.dags.project_name import project_name_of
from etna.etls.batches import AwaitBatches
from etna.etls.context import get_batch_range
from etna.etls.environment import etl_args
from etna.hooks.etna import EtnaHook, Metis, Folder, File
from etna.utils.inject import inject

# Increment this to force a reload of all metis cursor data.
metis_batch_loading_version = 1

def load_metis_files_and_folders_batch(
        metis: Metis,
        project_name: str, bucket_name: str,
        type: Union[Literal['file'], Literal['folder']]
    ) -> Union[List[File], List[Folder]]:
    context: Context = get_current_context()
    start, end = get_batch_range(context)

    result = metis.find(
        project_name, bucket_name, [
            dict(
                type=type,
                attribute='updated_at',
                predicate='>=',
                value=start.replace(tzinfo=timezone.utc).isoformat(timespec='seconds')
            ),
            dict(
                type=type,
                attribute='updated_at',
                predicate='<=',
                value=end.replace(tzinfo=timezone.utc).isoformat(timespec='seconds')
            ),
        ]
    )

    if type == 'file':
        return result.files

    return result.folders

#
# def load_metis_file_batches(project_name: str, bucket_name: str, start_date: datetime) -> DAG:
#     def load_batches_dag(batch_start_date, batch_end_date):
#         run_in_swarm('fetch_batch_from_metis', 'polyphemus_app', [
#             '/app/bin/polyphemus', 'etl', 'metis_file_tail_etl', 'find_batch', '--from-environment'
#         ], env=etl_args(dict(
#             limit=100,
#             project_bucket_pairs=[[project_name, bucket_name]],
#             cursor_env=dict(
#                 updated_at=batch_start_date,
#                 batch_end_at=batch_end_date,
#             )
#         )), output_json=True)
#
#     load_batches_dag.__name__ = f"metis_files_#{project_name}_#{bucket_name}_batch_loader"
#
#     return etl(
#         project_name,
#         start_date,
#         timedelta(minutes=5),
#         version=metis_batch_loading_version
#     )(load_batches_dag)
#
