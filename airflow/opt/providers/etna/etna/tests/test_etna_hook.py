import io
from datetime import datetime

import pytest
from airflow.models import Connection

from etna.hooks.etna import EtnaHook, File


@pytest.mark.vcr
def test_e2e_nonce_janus_to_metis(rsa_etna_connection: Connection):
    hook = EtnaHook(rsa_etna_connection.conn_id)
    with hook.metis("mvir1") as metis:
        list_result = metis.list_folder("mvir1", "data")
        assert list_result.files == []
        assert isinstance(list_result.folders[0].updated_at_datetime, datetime)

        list_result = metis.list_folder(
            "mvir1", "data", list_result.folders[0].folder_path
        )
        assert list_result.files != [] or list_result.folders != []


@pytest.mark.vcr
def test_e2e_token_to_metis(token_etna_connection: Connection):
    hook = EtnaHook(token_etna_connection.conn_id)
    with hook.metis("mvir1") as metis:
        list_result = metis.list_folder("mvir1", "data")
        assert list_result.files == []
        assert isinstance(list_result.folders[0].updated_at_datetime, datetime)

        list_result = metis.list_folder(
            "mvir1", "data", list_result.folders[0].folder_path
        )
        assert list_result.files != [] or list_result.folders != []

@pytest.mark.vcr
def test_e2e_upload_to_metis_and_download(token_etna_connection: Connection):
    hook = EtnaHook(token_etna_connection.conn_id)
    # Enough to trigger multiple upload batches
    payload = b'abcdef' * 2500000
    with hook.metis("ipi", read_only=False) as metis:
        # Adjust the file path to a unique file when re-recording.
        for upload in metis.upload_file('ipi', 'airflow_debug', 'test/a/file2', io.BytesIO(payload), len(payload)):
            print("Uploading...")
        print("Done!")

        with metis.open_file(upload.as_file, binary_mode=True) as download:
            result = download.read()

    assert payload == result