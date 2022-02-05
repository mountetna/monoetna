from datetime import datetime, timedelta
from unittest import mock

import pytest
from airflow import DAG
from airflow.decorators import task
from airflow.executors.debug_executor import DebugExecutor
from airflow.models import Connection, XCom
from airflow.utils.timezone import utc

from etna.dags import etl
from etna.etls.batches import AwaitBatches
from etna.etls.metis_files import load_metis_files_and_folders_batch
from etna.hooks.etna import EtnaHook
from etna.operators import run_in_container
from etna.tests.conftest import NotSoRandom


def run_dag(dag: DAG, execution_date: datetime, end_date: datetime):
    dag.run(
        executor=DebugExecutor(),
        start_date=execution_date,
        end_date=end_date,
        verbose=True,
    )

@pytest.mark.vcr
@mock.patch("tempfile._Random", NotSoRandom)
def test_metis_files_etl_e2e(token_etna_connection: Connection):
    hook = EtnaHook(token_etna_connection.conn_id)

    start = datetime(2022, 1, 1)

    def t(i=0):
        return (start + timedelta(minutes=i)).replace(tzinfo=utc)

    @etl('mvir1', datetime(2020, 1, 1), timedelta(minutes=5), 1)
    def load_metis_files_etl_dag():
        @task
        def load_files():
            with hook.metis() as metis:
                return load_metis_files_and_folders_batch(metis, 'data', 'folder')
        load_files()
    load_metis_files_etl_dag: DAG

    @etl('mvir1', datetime(2020, 1, 1), timedelta(minutes=5), 1)
    def process_some_metis_files_data():
        await_data_task = AwaitBatches(task_id='await_data', loader_dag=load_metis_files_etl_dag)

        output_task = run_in_container('output_batch', 'polyphemus_app', ['cat', '/inputs/batch'], output_json=True, docker_base_url="http://localhost:8085")
        output_task['/inputs/batch'] = await_data_task.output
    process_some_metis_files_data: DAG

    run_dag(load_metis_files_etl_dag, t(), t(4))
    files = XCom.get_one(dag_id=load_metis_files_etl_dag.dag_id, execution_date=t())
    assert len(files) > 0
    run_dag(load_metis_files_etl_dag, t(6), t(9))
    files = XCom.get_one(dag_id=load_metis_files_etl_dag.dag_id, execution_date=t(6))
    assert len(files) == 0

    run_dag(process_some_metis_files_data, t(0), t(4))
    files = XCom.get_one(dag_id=process_some_metis_files_data.dag_id, execution_date=t(0))
    assert len(files) > 0
