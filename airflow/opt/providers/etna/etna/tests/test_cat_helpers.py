import json
import re
from typing import Any
import pytz
from datetime import timedelta, datetime
from dateutil import parser
from unittest import mock
from pytest import raises

from airflow import DAG
from airflow.models.xcom_arg import XComArg
from airflow.models import Variable
from airflow.exceptions import AirflowException

from providers.etna.etna.etls.decorators import cat_etl
from providers.etna.etna.etls.cat import CatEtlHelpers, EtnaHook, CatHook, C4Hook
from providers.etna.etna.hooks.cat import SftpEntry, Cat
from providers.etna.etna.hooks.c4 import C4

from .test_metis_files_etl import run_dag, get_all_results

class MockSftpAttributes(object):
    def __init__(self, sftp_attrs: dict):
        self.sftp_attrs = sftp_attrs

    @property
    def filename(self):
        return self.sftp_attrs["filename"]

    @property
    def st_mtime(self):
        return self.sftp_attrs["st_mtime"]

    @property
    def st_size(self):
        return self.sftp_attrs["st_size"]

    @property
    def st_mode(self):
        return self.sftp_attrs["st_mode"]


def make_datetime(year, month, day):
    return datetime(year, month, day, 0, 0, 0, 0, pytz.UTC)


def timestamp(year, month, day):
    return int(make_datetime(year, month, day).timestamp())


def mock_tail():
    return [
        SftpEntry(MockSftpAttributes({
            "filename": "123.txt",
            "st_size": 1,
            "st_mtime": timestamp(2022, 1, 1)
        }), "parent/child/grandchild"),
        SftpEntry(MockSftpAttributes({
            "filename": "abc.exe",
            "st_size": 2,
            "st_mtime": timestamp(2022, 2, 1)
        }), ""),
        SftpEntry(MockSftpAttributes({
            "filename": "xyz.exe",
            "st_size":5,
            "st_mtime": timestamp(2022, 3, 1)
        }), "aunt/child/cousin")
    ]


def mock_var_get(system: str):
    if "c4" == system:
        return {
            "parent/child/grandchild/123.txt": f"1-{timestamp(2022, 1, 1)}",
            "something.txt": f"1-{timestamp(2022, 1, 1)}",
            "other_thing.txt": f"2-{timestamp(2022, 3, 1)}",
            "other_things.txt": f"3-{timestamp(2022, 5, 1)}",
            "folder/yet_another_thing.txt": f"5-{timestamp(2022, 6, 2)}",
            "folder/another_thing.txt": f"6-{timestamp(2022, 1, 2)}"
        }

    return {
        "abc.exe": f"2-{timestamp(2022, 2, 1)}",
        "other_thing.txt": f"0-{timestamp(2022, 3, 1)}",
        "other_things.txt": f"3-{timestamp(2022, 5, 1)}",
        "folder/yet_another_thing.txt": f"5-{timestamp(2022, 6, 2)}",
        "folder/another_thing.txt": f"6-{timestamp(2022, 1, 2)}"
    }


@mock.patch('providers.etna.etna.etls.decorators.load_cat_files_batch', side_effect=[mock_tail(), []])
def test_cat_files_etl_filter_file_name(mock_load, reset_db, ssh_cat_connection, ssh_c4_connection):
    record_matches_task_id: str = ""

    @cat_etl(version=1)
    def test_ingesting_cat_files_filter_file_name(helpers: CatEtlHelpers, tail_files: XComArg):
        nonlocal record_matches_task_id
        matches = helpers.filter_files(
            tail_files,
            file_name_regex=re.compile(r".*\.exe$")
        )

        record_matches_task_id = matches.operator.task_id

    test_ingesting_cat_files_filter_file_name: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_cat_files_filter_file_name, start_date, end_date)

    results = get_all_results(
        end_date, test_ingesting_cat_files_filter_file_name, record_matches_task_id
    )

    assert len(results) == 2
    assert results[0].name == "abc.exe"
    assert results[1].name == "xyz.exe"

@mock.patch('providers.etna.etna.etls.decorators.load_cat_files_batch', side_effect=[mock_tail(), []])
def test_cat_files_etl_filter_folder_path(mock_load, reset_db, ssh_cat_connection, ssh_c4_connection):
    record_matches_task_id: str = ""

    @cat_etl(version=1)
    def test_ingesting_cat_files_filter_folder_path(helpers: CatEtlHelpers, tail_files: XComArg):
        nonlocal record_matches_task_id
        matches = helpers.filter_files(
            tail_files,
            folder_path_regex=re.compile(r".*\/child\/.*")
        )

        record_matches_task_id = matches.operator.task_id

    test_ingesting_cat_files_filter_folder_path: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_cat_files_filter_folder_path, start_date, end_date)

    results = get_all_results(
        end_date, test_ingesting_cat_files_filter_folder_path, record_matches_task_id
    )

    assert len(results) == 2
    assert results[0].name == '123.txt'
    assert results[1].name == 'xyz.exe'


@mock.patch('providers.etna.etna.etls.decorators.load_cat_files_batch', side_effect=[mock_tail(), []])
def test_cat_files_etl_filter_file_name_and_folder_path(mock_load, reset_db, ssh_cat_connection, ssh_c4_connection):
    record_matches_task_id: str = ""

    @cat_etl(version=1)
    def test_ingesting_cat_files_filter_file_name_and_folder_path(helpers: CatEtlHelpers, tail_files: XComArg):
        nonlocal record_matches_task_id
        matches = helpers.filter_files(
            tail_files,
            file_name_regex=re.compile(r".*exe$"),
            folder_path_regex=re.compile(r".*\/child\/.*")
        )

        record_matches_task_id = matches.operator.task_id

    test_ingesting_cat_files_filter_file_name_and_folder_path: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_cat_files_filter_file_name_and_folder_path, start_date, end_date)

    results = get_all_results(
        end_date, test_ingesting_cat_files_filter_file_name_and_folder_path, record_matches_task_id
    )

    assert len(results) == 1
    assert results[0].name == 'xyz.exe'


def test_sftp_entry():
    entry_attrs_1 = MockSftpAttributes({
        "filename": "testfile.txt",
        "st_size": 42,
        "st_mtime": timestamp(2022, 1, 1),
        "st_mode": 33204 # is 0o100664, which is -rw-rw-r--
    })

    entry_attrs_2 = MockSftpAttributes({
        "filename": "A_FOLDER",
        "st_size": 1,
        "st_mtime": timestamp(2022, 2, 1),
        "st_mode": 16893 # is 0o40775, which is drwxrwxr-x
    })

    entry = SftpEntry(entry_attrs_1, "/root/child/leaf")
    assert entry.is_dir() == False
    assert entry.is_file() == True
    assert entry.size == 42
    assert entry.full_path == "/root/child/leaf/testfile.txt"
    assert entry.is_in_range() == True
    assert entry.is_in_range(make_datetime(2022, 3, 1)) == False
    assert entry.is_in_range(None, make_datetime(2021, 1, 1)) == False
    assert entry.is_in_range(make_datetime(2021, 1, 1)) == True
    assert entry.is_in_range(None, make_datetime(2022, 2, 1)) == True
    assert entry.is_in_range(make_datetime(2021, 1, 1), make_datetime(2022, 2, 1)) == True
    assert entry.is_in_range(make_datetime(2021, 1, 1), make_datetime(2021, 2, 1)) == False

    entry = SftpEntry(entry_attrs_2, "/root/child/leaf")
    assert entry.is_dir() == True
    assert entry.is_file() == False
    assert entry.size == 1
    assert entry.full_path == "/root/child/leaf/A_FOLDER"

mock_etna_hook = mock.Mock()
mock_metis = mock.Mock()
mock_c4 = mock.Mock()
mock_c4_hook = mock.Mock()

mock_retrieve_file = mock.Mock()

mock_cat_hook = mock.Mock()
mock_cat = mock.Mock()
mock_sftp = mock.Mock()
mock_c4_sftp = mock.Mock()
mock_listdir_attr = mock.Mock()
mock_cat_mark_ingested = mock.Mock()
mock_cat_update_cursor = mock.Mock()

def set_up_mocks():
    # Mock all this stuff so the context managers in ingest_to_metis
    #   are correctly handled
    mock_etna_hook.reset_mock()
    mock_metis.reset_mock()
    mock_c4.reset_mock()
    mock_c4_hook.reset_mock()
    mock_retrieve_file.reset_mock()
    mock_cat_hook.reset_mock()
    mock_cat.reset_mock()
    mock_sftp.reset_mock()
    mock_c4_sftp.reset_mock()
    mock_listdir_attr.reset_mock()
    mock_cat_mark_ingested.reset_mock()
    mock_cat_update_cursor.reset_mock()

    mock_c4_sftp.return_value.__enter__ = mock_c4_sftp
    mock_c4_sftp.return_value.__exit__ = mock_c4_sftp
    mock_c4.return_value.sftp = mock_c4_sftp
    mock_c4.return_value.upload_file.return_value = []
    mock_c4.return_value.__enter__ = mock_c4
    mock_c4.return_value.__exit__ = mock_c4
    mock_c4_hook.c4 = mock_c4

    mock_metis.return_value.upload_file.return_value = []
    mock_metis.return_value.__enter__ = mock_metis
    mock_metis.return_value.__exit__ = mock_metis
    mock_etna_hook.metis = mock_metis
    mock_listdir_attr.side_effect = [
        [
            MockSftpAttributes({
                "filename": "folder",
                "st_size": 1,
                "st_mtime": timestamp(2022, 1, 1),
                "st_mode": 16893 # is 0o40775, which is drwxrwxr-x
            }),
            MockSftpAttributes({
                "filename": "something.txt",
                "st_size": 1,
                "st_mtime": timestamp(2022, 2, 1),
                "st_mode": 33204 # is 0o100664, which is -rw-rw-r--
            }),
            MockSftpAttributes({
                "filename": "other_thing.txt",
                "st_size": 2,
                "st_mtime": timestamp(2022, 3, 1),
                "st_mode": 33204 # is 0o100664, which is -rw-rw-r--
            }),
            MockSftpAttributes({
                "filename": "other_things.txt",
                "st_size": 3,
                "st_mtime": timestamp(2022, 5, 1),
                "st_mode": 33204 # is 0o100664, which is -rw-rw-r--
            }),
        ], [
            MockSftpAttributes({
                "filename": "something_else.txt",
                "st_size": 1,
                "st_mtime": timestamp(2022, 3, 2),
                "st_mode": 33204 # is 0o100664, which is -rw-rw-r--
            }),
            MockSftpAttributes({
                "filename": "yet_another_thing.txt",
                "st_size": 1,
                "st_mtime": timestamp(2022, 6, 2),
                "st_mode": 33204 # is 0o100664, which is -rw-rw-r--
            }),
            MockSftpAttributes({
                "filename": "another_thing.txt",
                "st_size": 1,
                "st_mtime": timestamp(2022, 1, 2),
                "st_mode": 33204 # is 0o100664, which is -rw-rw-r--
            }),
        ]
    ]
    mock_sftp.return_value.listdir_attr = mock_listdir_attr
    mock_sftp.return_value.__enter__ = mock_sftp
    mock_sftp.return_value.__exit__ = mock_sftp
    mock_cat.return_value._root_path.return_value = 'ROOT'
    mock_cat.return_value.file_ingested_to_system.return_value = False
    mock_cat.return_value.sftp = mock_sftp

    mock_retrieve_file.return_value.__enter__ = mock_retrieve_file
    mock_retrieve_file.return_value.__exit__ = mock_retrieve_file
    mock_cat.return_value.mark_file_as_ingested = mock_cat_mark_ingested
    mock_cat.return_value.update_cursor = mock_cat_update_cursor
    mock_cat.return_value.retrieve_file = mock_retrieve_file
    mock_cat.return_value.__enter__ = mock_cat
    mock_cat.return_value.__exit__ = mock_cat
    mock_cat_hook.cat = mock_cat
    mock_cat_hook.connection.host = "host.development.local"

# Okay, these are not great tests, because I wound up mocking both the CAT
#   SFTP side as well as the C4 / Metis sides ... so it only tests the logic within
#   the helper, really. Probably fairly brittle.
@mock.patch('providers.etna.etna.etls.decorators.load_cat_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
def test_cat_files_etl_ingest_metis_remove_oligo_by_default(mock_etna, mock_load, reset_db):

    # Mock all this stuff so the context managers in ingest_to_metis
    #   are correctly handled
    set_up_mocks()

    @cat_etl(version=1, hook=mock_cat_hook, magic_string="12")
    def test_ingesting_cat_files_ingest_metis_remove_oligo_by_default(helpers: CatEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files, folder_path="foo")

    test_ingesting_cat_files_ingest_metis_remove_oligo_by_default: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_cat_files_ingest_metis_remove_oligo_by_default, start_date, end_date)

    mock_metis().upload_file.assert_any_call("triage", "waiting_room", "foo/parent/child/grandchild/3.txt", mock.ANY, 1)
    mock_retrieve_file.assert_called()


@mock.patch('providers.etna.etna.etls.decorators.load_cat_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(C4Hook, 'for_project', return_value=mock_c4_hook)
def test_cat_files_etl_ingest_c4_remove_oligo_by_default(mocked_c4, mock_load, reset_db):

    # Mock all this stuff so the context managers in ingest_to_metis
    #   are correctly handled
    set_up_mocks()

    @cat_etl(version=1, hook=mock_cat_hook, magic_string="12")
    def test_ingesting_cat_files_ingest_c4_remove_oligo_by_default(helpers: CatEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_c4(tail_files, folder_path="foo")

    test_ingesting_cat_files_ingest_c4_remove_oligo_by_default: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_cat_files_ingest_c4_remove_oligo_by_default, start_date, end_date)

    mock_c4().upload_file.assert_any_call(mock.ANY, "foo/parent/child/grandchild", "3.txt", mock.ANY, mock.ANY)
    mock_retrieve_file.assert_not_called()


@mock.patch('providers.etna.etna.etls.decorators.load_cat_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
def test_cat_files_etl_ingest_metis_do_not_remove_oligo(mock_etna, mock_load, reset_db):

    # Mock all this stuff so the context managers in ingest_to_metis
    #   are correctly handled
    set_up_mocks()

    @cat_etl(version=1, hook=mock_cat_hook, magic_string="12")
    def test_ingesting_cat_files_ingest_metis_do_not_remove_oligo(helpers: CatEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files, folder_path="foo", remove_magic_string=False)

    test_ingesting_cat_files_ingest_metis_do_not_remove_oligo: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_cat_files_ingest_metis_do_not_remove_oligo, start_date, end_date)

    mock_metis().upload_file.assert_any_call("triage", "waiting_room", "foo/parent/child/grandchild/123.txt", mock.ANY, 1)
    mock_retrieve_file.assert_called()


@mock.patch('providers.etna.etna.etls.decorators.load_cat_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(C4Hook, 'for_project', return_value=mock_c4_hook)
def test_cat_files_etl_ingest_c4_do_not_remove_oligo(mocked_c4, mock_load, reset_db):

    # Mock all this stuff so the context managers
    #   are correctly handled
    set_up_mocks()

    @cat_etl(version=1, hook=mock_cat_hook, magic_string="12")
    def test_ingesting_cat_files_ingest_c4_do_not_remove_oligo(helpers: CatEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_c4(tail_files, folder_path="foo", remove_magic_string=False)

    test_ingesting_cat_files_ingest_c4_do_not_remove_oligo: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_cat_files_ingest_c4_do_not_remove_oligo, start_date, end_date)

    mock_c4().upload_file.assert_any_call(mock.ANY, "foo/parent/child/grandchild", "123.txt", mock.ANY, mock.ANY)
    mock_retrieve_file.assert_not_called()


@mock.patch('providers.etna.etna.etls.decorators.load_cat_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
def test_cat_files_etl_ingest_metis_without_folder_path(mock_etna, mock_load, reset_db):

    set_up_mocks()

    @cat_etl(version=1, hook=mock_cat_hook)
    def test_ingesting_cat_files_ingest_metis_without_folder_path(helpers: CatEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files)

    test_ingesting_cat_files_ingest_metis_without_folder_path: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_cat_files_ingest_metis_without_folder_path, start_date, end_date)

    mock_metis().upload_file.assert_any_call("triage", "waiting_room", "host.development.local/parent/child/grandchild/123.txt", mock.ANY, 1)
    mock_retrieve_file.assert_called()


@mock.patch('providers.etna.etna.etls.decorators.load_cat_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(C4Hook, 'for_project', return_value=mock_c4_hook)
def test_cat_files_etl_ingest_c4_without_folder_path(mocked_c4, mock_load, reset_db):

    set_up_mocks()

    @cat_etl(version=1, hook=mock_cat_hook)
    def test_ingesting_cat_files_ingest_c4_without_folder_path(helpers: CatEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_c4(tail_files)

    test_ingesting_cat_files_ingest_c4_without_folder_path: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_cat_files_ingest_c4_without_folder_path, start_date, end_date)

    mock_c4().upload_file.assert_any_call(mock.ANY, "parent/child/grandchild", "123.txt", mock.ANY, mock.ANY)
    mock_retrieve_file.assert_not_called()


@mock.patch('providers.etna.etna.etls.decorators.load_cat_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
def test_cat_files_etl_ingest_metis_with_project_bucket(mock_etna, mock_load, reset_db):

    set_up_mocks()

    @cat_etl(version=1, hook=mock_cat_hook)
    def test_ingesting_cat_files_ingest_metis_with_project_bucket(helpers: CatEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files, project_name="test", bucket_name="bucket")

    test_ingesting_cat_files_ingest_metis_with_project_bucket: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_cat_files_ingest_metis_with_project_bucket, start_date, end_date)

    mock_metis().upload_file.assert_any_call("test", "bucket", "host.development.local/parent/child/grandchild/123.txt", mock.ANY, 1)
    mock_retrieve_file.assert_called()


@mock.patch('providers.etna.etna.etls.decorators.load_cat_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
@mock.patch.object(Variable, 'get', side_effect=[mock_var_get("c4"), mock_var_get("metis")])
@mock.patch.object(Variable, 'set')
def test_ingest_updates_cursor(mock_set, mock_get, mock_etna, mock_load, reset_db):

    set_up_mocks()

    @cat_etl(version=1, hook=mock_cat_hook)
    def test_ingest_updates_cursor(helpers: CatEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files)

    test_ingest_updates_cursor: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingest_updates_cursor, start_date, end_date)

    mock_cat_mark_ingested.assert_called()
    # 3 files from mock_tail should get ingested since one is in mock_var_get("metis")
    #   and a different one is in mock_var_get("c4") -- so no tail()
    #   file has been ingested into both systems.
    assert len(mock_cat_mark_ingested.call_args_list) == 3
    mock_cat_update_cursor.assert_called()


@mock.patch.object(Variable, 'get', side_effect=[mock_var_get("c4"), mock_var_get("metis")])
def test_tail(mock_var, reset_db):
    set_up_mocks()

    cat = Cat(mock_cat_hook)

    cat.sftp = mock_sftp
    cat._root_path = lambda: ""

    results = cat.tail(
        magic_string=re.compile(".*thing.*"),
        ignore_directories=['folder']
    )

    # Only two files are not in cursor or been updated, and not in an ignored directory
    #       but do match the magic string
    assert len(results) == 2
    assert [f.name for f in results].sort() == ['other_thing.txt', 'something.txt'].sort()


@mock.patch.object(Variable, 'get', side_effect=[mock_var_get("c4"), mock_var_get("metis")])
def test_mark_ingested_updates_cursor(mock_get, reset_db):
    set_up_mocks()

    cat = Cat(mock_cat_hook)

    cat.sftp = mock_sftp
    cat._root_path = lambda: ""

    test_file = SftpEntry(MockSftpAttributes({
                "filename": "another_thing.txt",
                "st_size": 1,
                "st_mtime": timestamp(2022, 1, 2),
                "st_mode": 33204 # is 0o100664, which is -rw-rw-r--
            }), "parent/child/grandchild")

    cat.cursors = {}
    cat.mark_file_as_ingested("c4", test_file)
    assert test_file.full_path in cat.cursors["c4"]
    assert test_file.hash == cat.cursors["c4"][test_file.full_path]


@mock.patch.object(Variable, 'get', side_effect=[mock_var_get("c4"), mock_var_get("metis")])
@mock.patch.object(Variable, 'set')
def test_update_cursor_saves_variable(mock_set, mock_get, reset_db):
    set_up_mocks()

    cat = Cat(mock_cat_hook)

    cat.sftp = mock_sftp
    cat._root_path = lambda: ""

    cat.cursors = {
        'c4': {
            '/file.txt': '1-2-3'
        },
        'metis': {}
    }
    cat.update_cursor("c4")
    mock_set.assert_called_with(
        'cat_ingest_cursor-c4',
        {
            '/file.txt': '1-2-3'
        },
        serialize_json=True)


def test_c4_validates_cat_config(reset_db):
    set_up_mocks()

    mock_c4_hook.connection.extra = json.dumps({
        "host_key": "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIOMqqnkVzrm0SdG6UOoqKLsabgH5C9okWi0dh2l9GKJl"
    })

    c4 = C4(mock_c4_hook)

    mock_cat_connection = mock.Mock()
    mock_cat_connection.host = "sub.domain.com"
    mock_cat_connection.login = "user123"
    mock_cat_connection.password = "pass_word"

    mock_cat_hook.connection = mock_cat_connection

    c4._validate_cat_config(mock_cat_hook)

    mock_cat_connection = mock.Mock()
    mock_cat_connection.host = "'sub.domain.com"
    mock_cat_connection.login = "user123#"
    mock_cat_connection.password = "pass_word'<>"

    mock_cat_hook.connection = mock_cat_connection

    with raises(AirflowException):
        c4._validate_cat_config(mock_cat_hook)

    mock_validate = mock.Mock()
    mock_ssh = mock.Mock()
    mock_stdout = mock.Mock()

    mock_stdout.channel.recv_exit_status = 0
    mock_ssh.exec_command = (mock_stdout, mock_stdout, mock_stdout)
    mock_ssh.return_value.__enter__ = mock_ssh
    mock_ssh.return_value.__exit__ = mock_ssh
    c4._validate_cat_config = mock_validate
    c4.ssh = mock_ssh

    test_file = SftpEntry(MockSftpAttributes({
                "filename": "another_thing.txt",
                "st_size": 1,
                "st_mtime": timestamp(2022, 1, 2),
                "st_mode": 33204 # is 0o100664, which is -rw-rw-r--
            }), "parent/child/grandchild")

    c4.upload_file(
        mock_c4_sftp,
        "foo",
        "file.txt",
        test_file,
        mock_cat_hook
    )

    mock_validate.assert_called()
