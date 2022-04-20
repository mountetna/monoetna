import re
from datetime import timedelta
from dateutil import parser
from unittest import mock

from airflow import DAG
from airflow.models.xcom_arg import XComArg


from providers.etna.etna.etls.decorators import box_etl
from providers.etna.etna.etls.box import BoxEtlHelpers
from providers.etna.etna.hooks.box import BoxFile
from providers.etna.etna.etls.box import EtnaHook

from .test_metis_files_etl import run_dag, get_all_results


def mock_tail():
    return [
        BoxFile(file_name="123.txt", folder_path="parent/child/grandchild"),
        BoxFile(file_name="abc.exe", folder_path=""),
        BoxFile(file_name="xyz.exe", folder_path="aunt/child/cousin")
    ]

@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
def test_metis_files_etl_filter_file_name(mock_load, reset_db, session_a):
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
    assert results[0].file_name == "abc.exe"
    assert results[1].file_name == "xyz.exe"

@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
def test_metis_files_etl_filter_folder_path(mock_load, reset_db, session_a):
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
    assert results[0].file_name == '123.txt'
    assert results[1].file_name == 'xyz.exe'


@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
def test_metis_files_etl_filter_file_name_and_folder_path(mock_load, reset_db, session_a):
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
    assert results[0].file_name == 'xyz.exe'


mock_etna_hook = mock.Mock()
mock_metis = mock.Mock()

mock_retrieve_file = mock.Mock()
mock_remove_file = mock.Mock()

# Okay, this is a really bad test, because I wound up mocking both the Box
#   FTP side as well as the Metis side ... so it only tests the logic within
#   the helper, really. Probably fairly brittle. Alternative might be to
#   record a VCR cassette, but not sure that works with FTPS?
# @mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
# @mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
# def test_metis_files_etl_ingest(mock_etna, mock_load, reset_db, session_a):

#     # Mock all this stuff so the context managers in ingest_to_metis
#     #   are correctly handled
#     mock_metis.return_value.upload_file.return_value = []
#     mock_metis.return_value.__enter__ = mock_metis
#     mock_metis.return_value.__exit__ = mock_metis
#     mock_etna_hook.metis = mock_metis
#     mock_box_hook = mock.Mock()
#     mock_box = mock.Mock()
#     mock_ftps = mock.Mock()
#     mock_ftps.return_value.__enter__ = mock_ftps
#     mock_ftps.return_value.__exit__ = mock_ftps
#     mock_box.return_value.ftps = mock_ftps

#     mock_retrieve_file.return_value.__enter__ = mock_retrieve_file
#     mock_retrieve_file.return_value.__exit__ = mock_retrieve_file
#     mock_box.return_value.file_size.return_value = 0
#     mock_box.return_value.retrieve_file = mock_retrieve_file
#     mock_box.return_value.remove_file = mock_remove_file
#     mock_box.return_value.__enter__ = mock_box
#     mock_box.return_value.__exit__ = mock_box
#     mock_box_hook.box = mock_box
#     mock_box_hook.connection.host = "host.development.local"

#     @box_etl("d_folder", version=1, hook=mock_box_hook)
#     def test_ingesting_metis_files_ingest(helpers: BoxEtlHelpers, tail_files: XComArg):
#         helpers.ingest_to_metis(tail_files, folder_path="foo", clean_up=True)

#     test_ingesting_metis_files_ingest: DAG

#     start_date = parser.parse("2022-01-01 00:00:00 +0000")
#     end_date = start_date + timedelta(days=1, minutes=1)
#     run_dag(test_ingesting_metis_files_ingest, start_date, end_date)

#     mock_metis().upload_file.assert_called()
#     mock_retrieve_file.assert_called()
#     mock_remove_file.assert_called()
