import re
import pytz
from datetime import timedelta, datetime
from dateutil import parser
from unittest import mock
from pytest import raises

from airflow import DAG
from airflow.models.xcom_arg import XComArg
from airflow.models import Variable


from providers.etna.etna.etls.decorators import box_etl
from providers.etna.etna.etls.box import BoxEtlHelpers
from providers.etna.etna.hooks.box import FtpEntry, Box
from providers.etna.etna.etls.box import EtnaHook

from .test_metis_files_etl import run_dag, get_all_results


def mock_tail():
    return [
        FtpEntry(("123.txt",{"size": 1, "modify": "20220101000000.000", "create": "20220301000000.000"}), "parent/child/grandchild"),
        FtpEntry(("abc.exe", {"size": 5, "modify": "20220101000000.000", "create": "20220301000000.000"}), ""),
        FtpEntry(("xyz.exe", {"size": 10, "modify": "20220101000000.000", "create": "20220301000000.000"}), "aunt/child/cousin")
    ]


def mock_var_get():
    return {
        "/something.txt": "1-20220201000000.000-20220101000000.000",
        "/other_thing.txt": "0-20220301000000.000-20220101000000.000",
        "/other_things.txt": '3-20220501000000.000-20220101000000.000',
        "/folder/yet_another_thing.txt": "5-20220602000000.000-20220101000000.000",
        "/folder/another_thing.txt": "6-20220102000000.000-20220101000000.000"
    }

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


def make_datetime(year, month, day):
    return datetime(year, month, day, 0, 0, 0, 0, pytz.UTC)


def test_ftp_entry():
    entry_dict_1 = {
        "modify": "20220101000000.000",
        "size": "42",
        "type": "file"
    }

    entry_dict_2 = {
        "modify": "20220201000000.000",
        "size": "42",
        "type": "dir"
    }

    entry = FtpEntry(("testfile.txt", entry_dict_1), "/root/child/leaf")
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

    entry = FtpEntry((".", entry_dict_2), "/root/child/leaf")
    assert entry.is_dir() == True
    assert entry.is_file() == False
    assert entry.is_dot() == True
    assert entry.size == 42
    assert entry.full_path == "/root/child/leaf/."

    entry = FtpEntry(("..", entry_dict_2), "/root/child/leaf")
    assert entry.is_dir() == True
    assert entry.is_file() == False
    assert entry.is_dot() == True
    assert entry.size == 42
    assert entry.full_path == "/root/child/leaf/.."

mock_etna_hook = mock.Mock()
mock_metis = mock.Mock()

mock_retrieve_file = mock.Mock()
mock_socket = mock.Mock()
mock_remove_file = mock.Mock()

mock_box_hook = mock.Mock()
mock_box = mock.Mock()
mock_ftps = mock.Mock()
mock_mlsd = mock.Mock()
mock_box_mark_ingested = mock.Mock()
mock_box_update_cursor = mock.Mock()

def set_up_mocks():
    # Mock all this stuff so the context managers in ingest_to_metis
    #   are correctly handled
    mock_etna_hook.reset_mock()
    mock_metis.reset_mock()
    mock_retrieve_file.reset_mock()
    mock_remove_file.reset_mock()
    mock_box_hook.reset_mock()
    mock_box.reset_mock()
    mock_ftps.reset_mock()
    mock_mlsd.reset_mock()
    mock_box_mark_ingested.reset_mock()
    mock_box_update_cursor.reset_mock()

    mock_metis.return_value.upload_file.return_value = []
    mock_metis.return_value.__enter__ = mock_metis
    mock_metis.return_value.__exit__ = mock_metis
    mock_etna_hook.metis = mock_metis
    mock_mlsd.side_effect = [
        [
            ('.', {'modify': '20220101000000.000', 'type': 'dir'}),
            ('folder', {'modify': '20220101000000.000', 'type': 'dir'}),
            ('something.txt', {'modify': '20220201000000.000', 'type': 'file', 'size': 1, 'create': '20220101000000.000'}),
            ('other_thing.txt', {'modify': '20220301000000.000', 'type': 'file', 'size': 2, 'create': '20220101000000.000'}),
            ('other_things.txt', {'modify': '20220501000000.000', 'type': 'file', 'size': 3, 'create': '20220101000000.000'})
        ], [
            ('.', {'modify': '20220101000000.000', 'type': 'dir'}),
            ('something_else.txt', {'modify': '20220302000000.000', 'type': 'file', 'size': 4, 'create': '20220101000000.000'}),
            ('yet_another_thing.txt', {'modify': '20220602000000.000', 'type': 'file', 'size': 5, 'create': '20220101000000.000'}),
            ('another_thing.txt', {'modify': '20220102000000.000', 'type': 'file', 'size': 6, 'create': '20220101000000.000'}),
        ]
    ]
    mock_ftps.return_value.mlsd = mock_mlsd
    mock_ftps.return_value.__enter__ = mock_ftps
    mock_ftps.return_value.__exit__ = mock_ftps
    mock_box.return_value.ftps = mock_ftps

    mock_retrieve_file.__enter__ = mock_retrieve_file
    mock_retrieve_file.__exit__ = mock_retrieve_file
    mock_socket.makefile.return_value = mock_retrieve_file
    mock_box.return_value.mark_file_as_ingested = mock_box_mark_ingested
    mock_box.return_value.update_cursor = mock_box_update_cursor
    mock_box.return_value.file_size.return_value = 0
    mock_box.return_value.retrieve_file.return_value = mock_socket
    mock_box.return_value.remove_file = mock_remove_file
    mock_box.return_value.__enter__ = mock_box
    mock_box.return_value.__exit__ = mock_box
    mock_box_hook.box = mock_box
    mock_box_hook.connection.host = "host.development.local"

# Okay, these are not great tests, because I wound up mocking both the Box
#   FTP side as well as the Metis side ... so it only tests the logic within
#   the helper, really. Probably fairly brittle. Alternative might be to
#   record a VCR cassette, but not sure that works with FTPS?
@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
def test_metis_files_etl_ingest_with_clean_up(mock_etna, mock_load, reset_db):

    set_up_mocks()

    @box_etl("d_folder", version=1, hook=mock_box_hook)
    def test_ingesting_metis_files_ingest_with_clean_up(helpers: BoxEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files, folder_path="foo", clean_up=True)

    test_ingesting_metis_files_ingest_with_clean_up: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_metis_files_ingest_with_clean_up, start_date, end_date)

    mock_metis().upload_file.assert_any_call("triage", "waiting_room", "foo/123.txt", mock.ANY, 1)
    mock_retrieve_file.assert_called()
    mock_remove_file.assert_called()


@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
def test_metis_files_etl_ingest_without_clean_up(mock_etna, mock_load, reset_db):

    # Mock all this stuff so the context managers in ingest_to_metis
    #   are correctly handled
    set_up_mocks()

    @box_etl("e_folder", version=1, hook=mock_box_hook)
    def test_ingesting_metis_files_ingest_without_clean_up(helpers: BoxEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files, folder_path="foo")

    test_ingesting_metis_files_ingest_without_clean_up: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_metis_files_ingest_without_clean_up, start_date, end_date)

    mock_metis().upload_file.assert_any_call("triage", "waiting_room", "foo/123.txt", mock.ANY, 1)
    mock_retrieve_file.assert_called()
    mock_remove_file.assert_not_called()


@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
def test_metis_files_etl_ingest_without_folder_path(mock_etna, mock_load, reset_db):

    set_up_mocks()

    @box_etl("f_folder", version=1, hook=mock_box_hook)
    def test_ingesting_metis_files_ingest_without_folder_path(helpers: BoxEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files)

    test_ingesting_metis_files_ingest_without_folder_path: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_metis_files_ingest_without_folder_path, start_date, end_date)

    mock_metis().upload_file.assert_any_call("triage", "waiting_room", "host.development.local/123.txt", mock.ANY, 1)
    mock_retrieve_file.assert_called()
    mock_remove_file.assert_not_called()


@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
def test_metis_files_etl_ingest_without_flatten(mock_etna, mock_load, reset_db):

    set_up_mocks()

    @box_etl("g_folder", version=1, hook=mock_box_hook)
    def test_ingesting_metis_files_ingest_without_flatten(helpers: BoxEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files, flatten=False)

    test_ingesting_metis_files_ingest_without_flatten: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_metis_files_ingest_without_flatten, start_date, end_date)

    mock_metis().upload_file.assert_any_call("triage", "waiting_room", "host.development.local/parent/child/grandchild/123.txt", mock.ANY, 1)
    mock_retrieve_file.assert_called()
    mock_remove_file.assert_not_called()


@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
def test_metis_files_etl_ingest_with_project_bucket(mock_etna, mock_load, reset_db):

    set_up_mocks()

    @box_etl("h_folder", version=1, hook=mock_box_hook)
    def test_ingesting_metis_files_ingest_with_project_bucket(helpers: BoxEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files, project_name="test", bucket_name="bucket")

    test_ingesting_metis_files_ingest_with_project_bucket: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_metis_files_ingest_with_project_bucket, start_date, end_date)

    mock_metis().upload_file.assert_any_call("test", "bucket", "host.development.local/123.txt", mock.ANY, 1)
    mock_retrieve_file.assert_called()
    mock_remove_file.assert_not_called()


@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
def test_metis_files_etl_ingest_with_split_folder_name(mock_etna, mock_load, reset_db):

    set_up_mocks()

    @box_etl("i_folder", version=1, hook=mock_box_hook)
    def test_ingesting_metis_files_ingest_with_split_folder_name(helpers: BoxEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files, flatten=False, split_folder_name="child")

    test_ingesting_metis_files_ingest_with_split_folder_name: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingesting_metis_files_ingest_with_split_folder_name, start_date, end_date)

    mock_metis().upload_file.assert_any_call("triage", "waiting_room", "host.development.local/grandchild/123.txt", mock.ANY, 1)
    mock_retrieve_file.assert_called()
    mock_remove_file.assert_not_called()


@mock.patch('providers.etna.etna.etls.decorators.load_box_files_batch', side_effect=[mock_tail(), []])
@mock.patch.object(EtnaHook, 'for_project', return_value=mock_etna_hook)
@mock.patch.object(Variable, 'get', side_effect=[mock_var_get()])
@mock.patch.object(Variable, 'set')
def test_ingest_updates_cursor(mock_set, mock_get, mock_etna, mock_load, reset_db):

    set_up_mocks()

    @box_etl("j_folder", version=1, hook=mock_box_hook)
    def test_ingest_updates_cursor(helpers: BoxEtlHelpers, tail_files: XComArg):
        helpers.ingest_to_metis(tail_files, flatten=False, split_folder_name="child")

    test_ingest_updates_cursor: DAG

    start_date = parser.parse("2022-01-01 00:00:00 +0000")
    end_date = start_date + timedelta(days=1, minutes=1)
    run_dag(test_ingest_updates_cursor, start_date, end_date)

    mock_box_mark_ingested.assert_called()
    # 3 files from mock_tail should all get ingested since none are in mock_var_get()
    assert len(mock_box_mark_ingested.call_args_list) == 3
    mock_box_update_cursor.assert_called()


@mock.patch.object(Variable, 'get', side_effect=[mock_var_get()])
def test_tail(mock_var, reset_db):
    set_up_mocks()

    box = Box(mock_box_hook)

    box.ftps = mock_ftps

    results = box.tail("ROOT")

    # Only two files are not in cursor or been updated
    assert len(results) == 2
    assert [f.name for f in results].sort() == ['other_thing.txt', 'something_else.txt'].sort()


def test_invalid_folder_name(reset_db):
    set_up_mocks()

    box = Box(mock_box_hook)

    with raises(ValueError):
        box.tail("/")


@mock.patch.object(Variable, 'get', side_effect=[mock_var_get()])
def test_mark_ingested_updates_cursor(mock_get, reset_db):
    set_up_mocks()

    box = Box(mock_box_hook)

    box.ftps = mock_ftps

    test_file = FtpEntry(("123.txt",{"size": 1, "modify": "20220101000000.000", "create": "20220301000000.000"}), "parent/child/grandchild")

    box.cursor = {}
    box.mark_file_as_ingested(test_file)
    assert test_file.full_path in box.cursor
    assert test_file.hash == box.cursor[test_file.full_path]


@mock.patch.object(Variable, 'get', side_effect=[mock_var_get()])
@mock.patch.object(Variable, 'set')
def test_update_cursor_saves_variable(mock_set, mock_get, reset_db):
    set_up_mocks()

    box = Box(mock_box_hook)

    box.ftps = mock_ftps

    box.cursor = {
        '/file.txt': '1-2-3'
    }
    box.update_cursor()
    mock_set.assert_called_with(
        'box_ingest_cursor',
        {
            '/file.txt': '1-2-3'
        },
        serialize_json=True)