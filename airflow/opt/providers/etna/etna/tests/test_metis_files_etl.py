import re
from datetime import datetime, timedelta
from typing import Literal, Union, List

from airflow import DAG
from airflow.executors.debug_executor import DebugExecutor
from airflow.models import Connection, XCom

from etna import metis_etl, MetisEtlHelpers
from etna.etls.etl_task_batching import LOWEST_BOUND

from etna.hooks.etna import Folder, File, EtnaHook
from etna.etls.metis import filter_by_record_directory


def run_dag(dag: DAG, execution_date: datetime, end_date: datetime):
    dag.run(
        executor=DebugExecutor(),
        start_date=execution_date,
        end_date=end_date,
        verbose=True,
    )


def get_all_results(end_date: datetime, dag: DAG, task_id: str):
    xcoms = XCom.get_many(
        execution_date=end_date,
        dag_ids=[dag.dag_id],
        task_ids=[task_id],
        include_prior_dates=True,
    )
    xcoms = xcoms.with_entities(XCom.value)
    return [y for x in xcoms for y in XCom.deserialize_value(x)]


def test_filter_by_root_directory():
    result = filter_by_record_directory(
        [
            Folder(folder_path="bulk_RNASeq/raw/abcdef"),
            File(file_path="bulk_RNASeq/raw/abcdef/myfile.txt"),
            File(file_path="bulk_RNASeq/raw/abcdef2/myfile.txt"),
        ],
        re.compile(r"bulk_RNASeq/raw/[^/]*"),
        "model",
    )

    assert result[0].root_path == "bulk_RNASeq/raw/abcdef"
    assert result[0].record_name == "abcdef"
    assert result[0].match_subpath == ""
    assert result[0].match_full_path == "bulk_RNASeq/raw/abcdef"
    assert result[0].folder_path == "bulk_RNASeq/raw/abcdef"
    assert result[0].match_file is None
    assert result[0].match_folder is not None

    assert result[1].root_path == "bulk_RNASeq/raw/abcdef2"
    assert result[1].record_name == "abcdef2"
    assert result[1].match_subpath == "myfile.txt"
    assert result[1].match_full_path == "bulk_RNASeq/raw/abcdef2/myfile.txt"
    assert result[1].folder_path == "bulk_RNASeq/raw/abcdef2"
    assert result[1].match_file is not None
    assert result[1].match_folder is None


def find_batch_start(
    hook: EtnaHook, project_name: str, bucket_name: str, folder_path: str
):
    with hook.metis(project_name) as metis:
        folders = metis.list_folder(
            project_name, bucket_name, folder_path=folder_path
        ).folders
        folders = sorted(folders, key=lambda f: f.updated_at_datetime)
        return folders[0].updated_at_datetime


# @mock.patch("tempfile._Random", NotSoRandom)
def test_metis_files_etl_e2e(token_etna_connection: Connection):
    hook = EtnaHook(token_etna_connection.conn_id)

    record_matches_task_id: str = ""

    @metis_etl("mvir1", "data", 1, hook=hook)
    def test_loading_metis_files(helpers: MetisEtlHelpers):
        nonlocal record_matches_task_id
        matches = helpers.find_record_folders(
            "rna_seq", re.compile(r"^bulk_RNASeq/raw/[^/]*")
        )
        matches = helpers.filter_by_timur(matches)
        matches = helpers.list_match_folders(matches)

        record_matches_task_id = matches.operator.task_id

    test_loading_metis_files: DAG

    start_date = find_batch_start(hook, "mvir1", "data", "bulk_RNASeq/raw")
    end_date = start_date + timedelta(days=1)
    run_dag(test_loading_metis_files, start_date, end_date)

    results = get_all_results(
        end_date, test_loading_metis_files, record_matches_task_id
    )
    assert len(results) > 0
