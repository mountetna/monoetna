import base64
import contextlib
import dataclasses
import hashlib
import io
import json
import logging
import os.path
import typing
import uuid
from datetime import datetime
from inspect import isgenerator
from typing import Dict, Optional, List
import re
from urllib.parse import quote

from airflow import DAG
from airflow.operators.python import get_current_context
from requests.adapters import HTTPAdapter
from serde import serialize, deserialize

import cached_property
import dateutil
import requests
from airflow.exceptions import AirflowException
from airflow.hooks.base import BaseHook
from airflow.models import Connection
from cryptography.hazmat.primitives import serialization, hashes
from cryptography.hazmat.primitives.asymmetric import padding
from requests import HTTPError, Session
from requests.auth import AuthBase
from serde.json import from_json, to_json
from urllib3 import Retry

from etna.dags.project_name import get_project_name
from etna.hooks.keys import prepared_key_from
from etna.utils.iterables import batch_iterable
from etna.utils.multipart import encode_as_multipart
from etna.utils.streaming import iterable_to_stream


class EtnaHook(BaseHook):
    """
    Etna Client Hook that can generate task tokens
    """

    conn_name_attr = "etna_conn_id"
    default_conn_name = "etna_default"
    conn_type = "etna"
    hook_name = "Etna Connection"
    root_dir = "/opt/airflow/dags/repos"

    remote_path: str
    local_path: str
    etna_conn_id: str
    key: Optional[str]

    # https://github.com/apache/airflow/blob/main/airflow/customized_form_field_behaviours.schema.json
    @staticmethod
    def get_ui_field_behaviour() -> Dict:
        return {
            "hidden_fields": ["port", "extra"],
            "relabeling": {
                "password": "RSA key or Janus Token",
                "host": "environment host postfix",
                "schema": "rsa or token",
                "login": "email of janus account (necessary for rsa)",
            },
            "placeholders": {
                "host": ".ucsf.edu",
                "schema": "token",
                "login": "you@ucsf.edu",
            },
        }

    @staticmethod
    def get_connection_form_widgets() -> Dict:
        """
        Returns dictionary of widgets to be added for the hook to handle extra values.

        If you have class hierarchy, usually the widgets needed by your class are already
        added by the base class, so there is no need to implement this method. It might
        actually result in warning in the logs if you try to add widgets that have already
        been added by the base class.

        Note that values of Dict should be of wtforms.Field type. It's not added here
        for the efficiency of imports.

        """
        return {}

    def __init__(self, etna_conn_id: str) -> None:
        super().__init__()
        self.etna_conn_id = etna_conn_id

    @classmethod
    def for_project(cls, project_name: Optional[str] = None):
        if project_name is None:
            project_name = get_project_name()
        return cls(f"etna_{project_name}")

    def get_conn(self) -> Connection:
        return self.get_connection(self.etna_conn_id)

    @cached_property.cached_property
    def connection(self) -> Connection:
        return self.get_conn()

    def get_hostname(self, host_prefix: str):
        return f"{host_prefix}{self.connection.host or '.ucsf.edu'}"

    def get_token_auth(self) -> AuthBase:
        schema = self.connection.schema or "token"
        if schema == "rsa":
            with prepared_key_from(self.connection, True) as key_file:
                with self.get_client(None) as session:
                    nonce = Janus(session, self.get_hostname("janus")).nonce()
                return SigAuth(key_file, self.connection.login, nonce)
        elif schema == "token":
            token = self.connection.password or ""
            return TokenAuth(token.encode("utf8"))
        else:
            raise AirflowException(
                f"Connection {self.connection.conn_id} has invalid schema: '{self.connection.schema}'"
            )

    def get_task_auth(
        self, project_name: Optional[str] = None, read_only=True
    ) -> "TokenAuth":
        if project_name is None:
            dag: DAG = get_current_context()["dag"]
            project_name = get_project_name(dag)
        token_auth = self.get_token_auth()
        with self.get_client(token_auth) as session:
            # Unfortunately, 'read only' here becomes a viewer token, which does not have
            # read permission to many buckets.  Unfortunately, read_only does not mean 'can read',
            # it means 'cannot write, but also likely cannot read.'
            # We are forced then to always create read_only=False tokens for now.
            token = Janus(session, self.get_hostname("janus")).generate_token(
                project_name, False
            )
        return TokenAuth(token, project_name)

    @contextlib.contextmanager
    def get_client(
        self, auth: Optional[AuthBase]
    ) -> typing.ContextManager[requests.Session]:
        yield EtnaSession(auth=auth, hooks={ response: [EtnaHook.assert_status] })

    @contextlib.contextmanager
    def metis(
        self, project_name: Optional[str] = None, read_only=True
    ) -> typing.ContextManager["Metis"]:
        auth = self.get_task_auth(project_name, read_only)
        with self.get_client(auth) as session:
            yield Metis(session, self.get_hostname("metis"))

    @contextlib.contextmanager
    def janus(
        self, project_name: Optional[str] = None, read_only=True
    ) -> typing.ContextManager["Janus"]:
        auth = self.get_task_auth(project_name, read_only)
        with self.get_client(auth) as session:
            yield Janus(session, self.get_hostname("janus"))

    @contextlib.contextmanager
    def magma(
        self, project_name: Optional[str] = None, read_only=True
    ) -> typing.ContextManager["Magma"]:
        auth = self.get_task_auth(project_name, read_only)
        with self.get_client(auth) as session:
            yield Magma(session, self.get_hostname("magma"))

    @staticmethod
    def assert_status(response: requests.Response, *args, **kwds):
        if 200 <= response.status_code < 300:
            return

        error_message = response.reason

        try:
            err_json = response.json()
            errors = err_json.get("errors", [])
            error = err_json.get("error")

            if error:
                errors.append(error)

            error_message = ", ".join(errors)
            error_message = (
                f"Request failed with {response.status_code}: {error_message}"
            )
        except (TypeError, ValueError, json.JSONDecodeError):
            pass

        raise HTTPError(error_message, response=response)

    def generate_task_token(self):
        pass



@serialize
@deserialize
@dataclasses.dataclass
class UploadResponse:
    current_byte_position: int = 0
    url: str = ""
    next_blob_size: int = 0

    @property
    def upload_path(self):
        return re.compile(r"^https://[^/]*?/").sub("/", self.url)


@dataclasses.dataclass
class Upload:
    file: typing.IO
    file_size: int
    upload_path: str

    cur: int = 0

    last_bytes: Optional[bytes] = None
    read_position: int = 0

    INITIAL_BLOB_SIZE = 2 ** 10
    MAX_BLOB_SIZE = 2 ** 22
    ZERO_HASH = "d41d8cd98f00b204e9800998ecf8427e"

    @property
    def as_magma_file_attribute(self) -> Optional[MagmaFileEntry]:
        metis_url = self.as_metis_url
        if metis_url:
            return MagmaFileEntry(
                path=metis_url, original_filename=os.path.basename(metis_url)
            )
        return None

    UPLOAD_PATH_REGEX = re.compile(r"/([^/]+)/upload/([^/]+)/?([^?]*)")

    @property
    def as_metis_url(self) -> Optional[str]:
        match = self.UPLOAD_PATH_REGEX.match(self.upload_path)
        if match:
            project_name = match.group(1)
            bucket_name = match.group(2)
            path = match.group(3)
            return f"metis://{project_name}/{bucket_name}/{path}"
        return None

    @property
    def as_file(self) -> Optional[File]:
        match = self.UPLOAD_PATH_REGEX.match(self.upload_path)
        if match:
            project_name = match.group(1)
            bucket_name = match.group(2)
            path = match.group(3)
            return File(file_name=os.path.basename(path), file_path=path, project_name=project_name, bucket_name=bucket_name)
        return None

    @property
    def next_blob_size(self):
        return min(
            self.MAX_BLOB_SIZE if self.cur > 0 else self.INITIAL_BLOB_SIZE,
            self.file_size - self.cur,
        )

    @property
    def is_complete(self):
        return self.cur >= self.file_size

    @property
    def next_blob_hash(self):
        b = self.read_next_blob()
        if not b:
            return self.ZERO_HASH

        return hashlib.md5(b).hexdigest()

    def resume_from(self, upload_response: UploadResponse):
        self.cur = upload_response.current_byte_position

    def advance_position(self):
        self.cur += self.next_blob_size

    def read_next_blob(self) -> bytes:
        next_left = self.cur
        next_right = self.cur + self.next_blob_size

        # Advance the position, if we've resumed to a specific place that needs skipping to.
        # We are merely reading data to catch up to our expected left position.
        if next_right < self.read_position:
            if self.file.seekable():
                self.file.seek(next_left)
                self.read_position = next_left
            else:
                raise ValueError(
                    "Upload source stream requires restart, upload has failed."
                )
        elif self.read_position < next_left:
            bytes_to_read = next_left - self.read_position
            data = self.file.read(bytes_to_read)
            self.read_position += len(data)
            if not len(data) < bytes_to_read:
                raise ValueError("Unexpected EOF while reading source stream.")

        if next_right == self.read_position:
            return self.last_bytes

        if next_left != self.read_position:
            if self.file.seekable():
                self.file.seek(self.cur)
                self.read_position = self.cur
            else:
                raise ValueError(
                    "Alignment error, upload needs to restart but source is not seekable."
                )

        data: List[bytes] = []
        while self.read_position < next_right:
            expected_bytes = next_right - self.read_position
            chunk = self.file.read(expected_bytes)
            self.read_position += len(chunk)
            if len(chunk) < expected_bytes:
                raise ValueError("Unexpected EOF while reading source stream")
            data.append(chunk)

        self.last_bytes = b"".join(data)
        return self.last_bytes


@serialize
@deserialize
@dataclasses.dataclass
class UploadStartRequest:
    file_size: int = 0
    action: str = "start"
    metis_uid: str = ""
    next_blob_size: int = 0
    upload_path: str = ""
    next_blob_hash: str = ""
    reset: bool = False


@dataclasses.dataclass
class UploadBlobRequest:
    file_size: int = 0
    action: str = "blob"
    metis_uid: str = ""
    blob_data: bytes = b""
    upload_path: str = ""
    next_blob_size: int = 0
    next_blob_hash: str = ""
    current_byte_position: int = 0

    def to_dict(self):
        return self.__dict__


