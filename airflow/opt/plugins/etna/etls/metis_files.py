import functools
from datetime import timedelta, datetime
from typing import Union

from airflow import DAG

from etna import etl, run_in_swarm
from etna.dags.project_name import project_name_of
from etna.etls.batches import AwaitBatches
from etna.etls.environment import etl_args
from etna.utils.inject import inject

# Increment this to force a reload of all metis cursor data.
metis_batch_loading_version = 1

def process_metis_files(loading_dag: DAG, start_date: datetime, interval: timedelta, version: Union[int, str]):
    version = f"{version}_{metis_batch_loading_version}"

    def instantiate_dag(fn):
        @functools.wraps(fn)
        def load_batch_first():
            return inject(fn, dict(load_batch=AwaitBatches(loading_dag, 'updated_at', task_id='load_metis_file_batch')))

        return etl(
            project_name=project_name_of(loading_dag),
            start_date=start_date,
            interval=interval,
            version=version,
        )(load_batch_first)

    return instantiate_dag

def load_metis_file_batches(project_name: str, bucket_name: str, start_date: datetime) -> DAG:
    def load_batches_dag(batch_start_date, batch_end_date):
        run_in_swarm('fetch_batch_from_metis', 'polyphemus_app', [
            '/app/bin/polyphemus', 'etl', 'metis_file_tail_etl', 'find_batch', '--from-environment'
        ], env=etl_args(dict(
            limit=100,
            project_bucket_pairs=[[project_name, bucket_name]],
            cursor_env=dict(
                updated_at=batch_start_date,
                batch_end_at=batch_end_date,
            )
        )), output_json=True)

    load_batches_dag.__name__ = f"metis_files_#{project_name}_#{bucket_name}_batch_loader"

    return etl(
        project_name,
        start_date,
        timedelta(minutes=5),
        version=metis_batch_loading_version
    )(load_batches_dag)

