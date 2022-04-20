import re
from datetime import timedelta
from dateutil import parser
from unittest import mock

from airflow import DAG
from airflow.models.xcom_arg import XComArg


from providers.etna.etna.etls.decorators import box_etl
from providers.etna.etna.etls.box import BoxEtlHelpers
from providers.etna.etna.hooks.box import FtpEntry
from providers.etna.etna.etls.box import EtnaHook

from .test_metis_files_etl import run_dag, get_all_results


def mock_tail():
    return [
        FtpEntry(("123.txt",{}), "parent/child/grandchild"),
        FtpEntry(("abc.exe", {}), ""),
        FtpEntry(("xyz.exe", {}), "aunt/child/cousin")
    ]

@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
def test_metis_files_etl_filter_file_name(mock_load, reset_db):
    record_matches_task_id: str = ""

    @box_etl("a_folder", version=1)
    def test_ingesting_metis_files_filter_file_name(helpers: BoxEtlHelpers, tail_files: XComArg):
        nonlocal record_matches_task_id
        matches = helpers.filter_files(
            tail_files,
            file_name_regex=re.compile(r".*\.exe$")
        )

        record_matches_task_id = matches.operator.task_id

    test_ingesting_metis_files_filter_file_name: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_metis_files_filter_file_name, start_date, end_date)

    results = get_all_results(
        end_date, test_ingesting_metis_files_filter_file_name, record_matches_task_id
    )

    assert len(results) == 2
    assert results[0].name == "abc.exe"
    assert results[1].name == "xyz.exe"

@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
def test_metis_files_etl_filter_folder_path(mock_load, reset_db):
    record_matches_task_id: str = ""

    @box_etl("b_folder", version=1)
    def test_ingesting_metis_files_filter_folder_path(helpers: BoxEtlHelpers, tail_files: XComArg):
        nonlocal record_matches_task_id
        matches = helpers.filter_files(
            tail_files,
            folder_path_regex=re.compile(r".*\/child\/.*")
        )

        record_matches_task_id = matches.operator.task_id

    test_ingesting_metis_files_filter_folder_path: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_metis_files_filter_folder_path, start_date, end_date)

    results = get_all_results(
        end_date, test_ingesting_metis_files_filter_folder_path, record_matches_task_id
    )

    assert len(results) == 2
    assert results[0].name == '123.txt'
    assert results[1].name == 'xyz.exe'


@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
def test_metis_files_etl_filter_file_name_and_folder_path(mock_load, reset_db):
    record_matches_task_id: str = ""

    @box_etl("c_folder", version=1)
    def test_ingesting_metis_files_filter_file_name_and_folder_path(helpers: BoxEtlHelpers, tail_files: XComArg):
        nonlocal record_matches_task_id
        matches = helpers.filter_files(
            tail_files,
            file_name_regex=re.compile(r".*exe$"),
            folder_path_regex=re.compile(r".*\/child\/.*")
        )

        record_matches_task_id = matches.operator.task_id

    test_ingesting_metis_files_filter_file_name_and_folder_path: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_metis_files_filter_file_name_and_folder_path, start_date, end_date)

    results = get_all_results(
        end_date, test_ingesting_metis_files_filter_file_name_and_folder_path, record_matches_task_id
    )

    assert len(results) == 1
    assert results[0].name == 'xyz.exe'
