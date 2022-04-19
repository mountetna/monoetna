import re
from datetime import timedelta
from dateutil import parser
from typing import List
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

# @mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
# def test_metis_files_etl_filter_file_name(reset_db):
#     record_matches_task_id: str = ""

#     @box_etl("a_folder", version=1)
#     def test_ingesting_metis_files_filter_file_name(helpers: BoxEtlHelpers, tail_files: XComArg):
#         nonlocal record_matches_task_id
#         matches = helpers.filter_files(
#             tail_files,
#             file_name_regex=re.compile(r".*\.exe$")
#         )

#         record_matches_task_id = matches.operator.task_id

#     test_ingesting_metis_files_filter_file_name: DAG

#     start_date = parser.parse("2022-01-01 00:00:00 +0000")
#     end_date = start_date + timedelta(days=1, minutes=1)
#     run_dag(test_ingesting_metis_files_filter_file_name, start_date, end_date)

#     results = get_all_results(
#         end_date, test_ingesting_metis_files_filter_file_name, record_matches_task_id
#     )

#     assert len(results) == 1
#     assert results[0].file_name == "abc.exe"

# @mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
# def test_metis_files_etl_filter_folder_path(reset_db):
#     record_matches_task_id: str = ""

#     @box_etl("b_folder", version=1)
#     def test_ingesting_metis_files_filter_folder_path(helpers: BoxEtlHelpers, tail_files: XComArg):
#         nonlocal record_matches_task_id
#         matches = helpers.filter_files(
#             tail_files,
#             folder_path_regex=re.compile(r".*\/child\/.*")
#         )

#         record_matches_task_id = matches.operator.task_id

#     test_ingesting_metis_files_filter_folder_path: DAG

#     start_date = parser.parse("2022-01-01 00:00:00 +0000")
#     end_date = start_date + timedelta(days=1, minutes=1)
#     run_dag(test_ingesting_metis_files_filter_folder_path, start_date, end_date)

#     results = get_all_results(
#         end_date, test_ingesting_metis_files_filter_folder_path, record_matches_task_id
#     )

#     assert len(results) == 1
#     assert results[0].file_name == '123.txt'


# @mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
# def test_metis_files_etl_filter_file_name_and_folder_path(reset_db):
#     record_matches_task_id: str = ""

#     @box_etl("c_folder", version=1)
#     def test_ingesting_metis_files_filter_file_name_and_folder_path(helpers: BoxEtlHelpers, tail_files: XComArg):
#         nonlocal record_matches_task_id
#         matches = helpers.filter_files(
#             tail_files,
#             file_name_regex=re.compile(r".*exe$"),
#             folder_path_regex=re.compile(r".*\/child\/.*")
#         )

#         record_matches_task_id = matches.operator.task_id

#     test_ingesting_metis_files_filter_file_name_and_folder_path: DAG

#     start_date = parser.parse("2022-01-01 00:00:00 +0000")
#     end_date = start_date + timedelta(days=1, minutes=1)
#     run_dag(test_ingesting_metis_files_filter_file_name_and_folder_path, start_date, end_date)

#     results = get_all_results(
#         end_date, test_ingesting_metis_files_filter_file_name_and_folder_path, record_matches_task_id
#     )

#     assert len(results) == 1
#     assert results[0].file_name == 'xyz.exe'


mock_etna_hook = mock.Mock()
mock_metis = mock.Mock()
mock_etna_hook.metis.return_value = mock_metis

mock_retrieve_file = mock.Mock()
mock_remove_file = mock.Mock()

@mock.patch('providers.etna.etna.etls.box._retrieve_file', return_value=mock_retrieve_file)
@mock.patch('providers.etna.etna.etls.box._remove_file', return_value=mock_remove_file)
@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
def test_metis_files_etl_ingest(mock_etna, mock_load, mock_remove, mock_retrieve, reset_db):

    @box_etl("d_folder", version=1)
    def test_ingesting_metis_files_ingest(helpers: BoxEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files)

    test_ingesting_metis_files_ingest: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_metis_files_ingest, start_date, end_date)

    mock_metis.upload_file.assert_called()
    mock_retrieve.assert_called()
    mock_remove.assert_called()
