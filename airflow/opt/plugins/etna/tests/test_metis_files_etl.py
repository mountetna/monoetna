import re
from datetime import datetime, timedelta
from typing import List, Tuple
from unittest import mock

import pytest
from airflow import DAG
from airflow.decorators import task
from airflow.executors.debug_executor import DebugExecutor
from airflow.models import Connection, XCom
from airflow.utils.timezone import utc

from etna import etl
from etna.dags.decorators import metis_folder_etl
from etna.etls.batches import AwaitBatches
from etna.etls.metis import filter_by_root_directory, MatchedAtRoot, filter_by_exists_in_timur, link, \
    list_contents_of_matches
from etna.hooks.etna import EtnaHook, Folder, FoldersAndFilesResponse, File
from etna.operators import run_in_container
from etna.tests.conftest import NotSoRandom
from etna.xcom.etna_xcom import pickled


def run_dag(dag: DAG, execution_date: datetime, end_date: datetime):
    dag.run(
        executor=DebugExecutor(),
        start_date=execution_date,
        end_date=end_date,
        verbose=True,
    )

def test_filter_by_root_directory():
    result = filter_by_root_directory([
        Folder(folder_path="bulk_RNASeq/raw/abcdef"),
    ], re.compile(r'bulk_RNASeq/[^/]*'))

    assert result[0].root_path == 'bulk_RNASeq/raw'
    assert result[0].root_base_name == 'raw'
    assert result[0].match_subpath == 'abcdef'
    assert result[0].match_full_path == 'bulk_RNASeq/raw/abcdef'

@pytest.mark.vcr
@mock.patch("tempfile._Random", NotSoRandom)
def test_metis_files_etl_e2e(token_etna_connection: Connection):
    hook = EtnaHook(token_etna_connection.conn_id)

    start = datetime(2022, 1, 1)

    def t(i=0):
        return (start + timedelta(minutes=i)).replace(tzinfo=utc)

    @metis_folder_etl('mvir1', 'data', datetime(2020, 1, 1), 1, hook=hook)
    def load_metis_files_etl_dag(tail_folders, project_name, bucket_name):
        @task
        def find_rna_seq_record_folders(folders: List[Folder]):
            matching_record_folders = filter_by_root_directory(folders, re.compile(r'^bulk_RNASeq/raw/[^/]*'))
            with hook.magma() as magma:
                matching_record_folders = filter_by_exists_in_timur(magma, matching_record_folders, 'rna_seq')
            with hook.metis() as metis:
                return list_contents_of_matches(metis, matching_record_folders)

        @link('rna_seq', 'raw_fastq_files', dry_run=True)
        def link_raw_fastq_files(matches: List[Tuple[MatchedAtRoot, File]]):
            for match, file in matches:
                if re.compile(r'/.*\.fastq\.gz$').match(file.file_path):
                    yield match.root_base_name, file

        rna_seq_record_folders = find_rna_seq_record_folders(tail_folders)
        link_raw_fastq_files(rna_seq_record_folders)


    load_metis_files_etl_dag: DAG

    @etl('mvir1', datetime(2020, 1, 1), timedelta(minutes=5), 1)
    def process_with_docker():
        await_data_task = AwaitBatches(task_id='await_data', loader_dag_or_id=load_metis_files_etl_dag, loader_task_or_id='tail_folders')
        output_task = run_in_container('output_batch', 'polyphemus_app', ['cat', '/inputs/batch'], output_json=True, docker_base_url="http://localhost:8085")
        output_task['/inputs/batch'] = await_data_task.output
    process_with_docker: DAG

    run_dag(load_metis_files_etl_dag, t(), t(4))
    folders = XCom.get_one(dag_id=load_metis_files_etl_dag.dag_id, task_id='tail_folders', execution_date=t())
    assert len(folders) > 0
    run_dag(load_metis_files_etl_dag, t(6), t(9))
    folders = XCom.get_one(dag_id=load_metis_files_etl_dag.dag_id, task_id='tail_folders', execution_date=t(6))
    assert len(folders) == 0

    matches: List[MatchedAtRoot] = XCom.get_one(dag_id=load_metis_files_etl_dag.dag_id, task_id='find_record_folders', execution_date=t())
    assert len(matches) > 0
    assert matches[0].match_folder.folder_path.startswith('bulk_RNASeq/')
    assert matches[0].root_path == 'bulk_RNASeq/raw'
    assert matches[0].match_subpath == 'MYVIR-HS15-D0PBMC1-RSQ1'

    run_dag(process_with_docker, t(0), t(4))
    folders = XCom.get_one(dag_id=process_with_docker.dag_id, task_id='output_batch', execution_date=t(0))
    assert len(folders) > 0