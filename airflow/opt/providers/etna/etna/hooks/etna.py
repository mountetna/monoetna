import base64
import contextlib
import dataclasses
import hashlib
import io
import json
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
            token = Janus(session, self.get_hostname("janus")).generate_token(
                project_name, read_only
            )
        return TokenAuth(token, project_name)

    @contextlib.contextmanager
    def get_client(
        self, auth: Optional[AuthBase]
    ) -> typing.ContextManager[requests.Session]:
        with requests.Session() as session:
            session.mount(
                "https://",
                HTTPAdapter(
                    max_retries=Retry(
                        total=5,
                        backoff_factor=1,
                        connect=5,
                        read=3,
                        status_forcelist=[502, 503, 504],
                    )
                ),
            )

            if auth:
                session.auth = auth
            session.hooks["response"].append(EtnaHook.assert_status)
            yield session

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


class TokenAuth(AuthBase):
    token: bytes
    project_scope: Optional[str]

    def __init__(self, token: bytes, project_scope: Optional[str] = None):
        self.token = token
        self.project_scope = project_scope

    def __call__(self, r: requests.Request) -> requests.Request:
        r.headers["Authorization"] = f'Etna {self.token.decode("ascii")}'
        return r


class SigAuth(AuthBase):
    nonce: bytes
    email: str
    private_key: bytes

    def __init__(self, private_key_file: str, email: str, nonce: bytes):
        self.private_key = open(private_key_file, "rb").read()
        self.nonce = nonce
        self.email = email

    def __call__(self, r: requests.Request) -> requests.Request:
        key = serialization.load_pem_private_key(self.private_key, None)
        txt_to_sign: bytes = b".".join(
            [self.nonce, base64.b64encode(self.email.encode("ascii"))]
        )
        sig = base64.b64encode(
            key.sign(
                txt_to_sign,
                padding.PKCS1v15(
                    # If in the future the server changes to PSS
                    # mgf=padding.MGF1(hashes.SHA256()),
                    # salt_length=padding.PSS.MAX_LENGTH
                ),
                hashes.SHA256(),
            )
        )
        together: bytes = b".".join([txt_to_sign, sig])
        r.headers["Authorization"] = f'Signed-Nonce {together.decode("ascii")}'
        return r


def encode_path(*segments: str) -> str:
    return "/".join(quote(s) for s in segments)


class EtnaClientBase:
    session: requests.Session
    hostname: str

    def __init__(self, session: requests.Session, hostname: str):
        self.session = session
        self.hostname = hostname

    @property
    def auth_token(self) -> Optional[bytes]:
        if isinstance(self.session.auth, TokenAuth):
            return self.session.auth.token

    def prepare_url(self, *path: str):
        return f"https://{self.hostname}/{encode_path(*path)}"

    def prepare_url_unsafe(self, path: str):
        return f"https://{self.hostname}/{path}"

    def get_project_scope(self) -> Optional[str]:
        if isinstance(self.session.auth, TokenAuth):
            return self.session.auth.project_scope
        return None


class Janus(EtnaClientBase):
    def generate_token(
        self, project_name: str, ready_only=True, token_type="task"
    ) -> bytes:
        response = self.session.post(
            self.prepare_url("api", "tokens", "generate"),
            json=dict(
                token_type=token_type, project_name=project_name, ready_only=ready_only
            ),
        )
        return response.content

    def nonce(self) -> bytes:
        response = self.session.get(self.prepare_url("api", "tokens", "nonce"))
        return response.content


@serialize
@deserialize
@dataclasses.dataclass
class MagmaFileEntry:
    path: str = ""
    original_filename: str = ""


@serialize
@deserialize
@dataclasses.dataclass
class File:
    file_name: str = ""
    file_hash: str = ""
    updated_at: str = ""
    file_path: Optional[str] = ""
    folder_id: Optional[int] = None
    project_name: str = ""
    bucket_name: str = ""
    download_url: Optional[str] = ""

    @property
    def as_magma_file_attribute(self) -> MagmaFileEntry:
        return MagmaFileEntry(
            path=self.as_metis_url,
            original_filename=self.file_name or os.path.basename(self.file_path),
        )

    @property
    def as_metis_url(self):
        return f"metis://{self.project_name}/{self.bucket_name}/{self.file_path}"

    @property
    def updated_at_datetime(self) -> datetime:
        return dateutil.parser.isoparse(self.updated_at)


@serialize
@deserialize
@dataclasses.dataclass
class Folder:
    id: int = 0
    folder_path: str = ""
    project_name: str = ""
    bucket_name: str = ""
    folder_name: str = ""
    updated_at: str = ""

    @property
    def updated_at_datetime(self) -> datetime:
        return dateutil.parser.isoparse(self.updated_at)


@serialize
@deserialize
@dataclasses.dataclass
class FoldersAndFilesResponse:
    folders: List[Folder]
    files: List[File]

    def empty(self):
        return not self.folders and not self.files

    def extend(self, other: "FoldersAndFilesResponse"):
        self.folders.extend(other.folders)
        self.files.extend(other.files)


@serialize
@deserialize
@dataclasses.dataclass
class FoldersResponse:
    folders: List[Folder]


class Metis(EtnaClientBase):
    def touch_folder(
        self, project_name: str, bucket_name: str, folder_path: str
    ) -> FoldersResponse:
        response = self.session.get(
            self.prepare_url(project_name, "folder", "touch", bucket_name, folder_path)
        )
        return from_json(FoldersResponse, response.content)

    def list_folder(
        self, project_name: str, bucket_name: str, folder_path: Optional[str] = None
    ) -> FoldersAndFilesResponse:
        if folder_path:
            response = self.session.get(
                self.prepare_url(project_name, "list", bucket_name, folder_path)
            )
        else:
            response = self.session.get(
                self.prepare_url(project_name, "list", bucket_name)
            )
        return from_json(FoldersAndFilesResponse, response.content)

    def open_file(self, file: File) -> typing.ContextManager[io.BufferedReader]:
        return self.open_url(file.download_url)

    @contextlib.contextmanager
    def open_url(self, url: str) -> typing.ContextManager[io.BufferedReader]:
        response = self.session.get(url, stream=True)
        yield iterable_to_stream(response.iter_content(io.DEFAULT_BUFFER_SIZE))
        response.close()

    def authorize_upload(
        self, project_name: str, bucket_name: str, file_path: str
    ) -> "UploadResponse":
        response = self.session.post(
            self.prepare_url("authorize", "upload"),
            json=dict(
                project_name=project_name, bucket_name=bucket_name, file_path=file_path
            ),
        )
        response_obj = from_json(UploadResponse, response.content)
        return response_obj

    def start_upload(self, request: "UploadStartRequest") -> "UploadResponse":
        response = self.session.post(
            self.prepare_url_unsafe(request.upload_path),
            data=to_json(request),
            headers={"Content-Type": "application/json"},
        )
        response_obj = from_json(UploadResponse, response.content)
        return response_obj

    def upload_blob(self, request: "UploadBlobRequest"):
        response = self.session.post(
            self.prepare_url_unsafe(request.upload_path),
            files=encode_as_multipart(request),
        )
        response_obj = from_json(UploadResponse, response.content)
        return response_obj

    def upload_file(
        self,
        project_name: str,
        bucket_name: str,
        dest_path: str,
        file: typing.IO,
        size: Optional[int] = None,
        max_retries=3,
    ) -> typing.Iterable["Upload"]:
        if size is None:
            if file.seekable():
                file.seek(0, 2)
                size = file.tell()
                file.seek(0)
            else:
                raise ValueError(
                    "size was not specified, but file is streaming.  You must specify a size hint when using a stream."
                )

        if len(dest_path) > 1:
            self.create_folder(project_name, bucket_name, os.path.dirname(dest_path))
        authorization = self.authorize_upload(project_name, bucket_name, dest_path)
        upload = Upload(file=file, size=size, upload_path=authorization.upload_path)
        return self.upload_parts(upload, uuid.uuid4().hex, max_retries)

    def upload_parts(
        self, upload: "Upload", metis_uid: str, remaining_attempts: int, reset=False
    ) -> typing.Iterable["Upload"]:
        unsent_zero_byte_file = upload.cur == 0
        upload.resume_from(
            self.start_upload(
                UploadStartRequest(
                    upload_path=upload.upload_path,
                    file_size=upload.size,
                    next_blob_size=upload.next_blob_size,
                    next_blob_hash=upload.next_blob_hash,
                    metis_uid=metis_uid,
                    reset=reset,
                )
            )
        )

        while unsent_zero_byte_file and not upload.is_complete:
            try:
                self.upload_blob(
                    UploadBlobRequest(
                        upload_path=upload.upload_path,
                        next_blob_size=upload.next_blob_size,
                        next_blob_hash=upload.next_blob_hash,
                        blob_data=upload.read_next_blob(),
                        metis_uid=metis_uid,
                        current_byte_position=upload.cur,
                    )
                )

                unsent_zero_byte_file = False
            except HTTPError as e:
                if remaining_attempts > 1:
                    if e.response.status_code == 422:
                        yield from self.upload_parts(
                            upload, metis_uid, remaining_attempts - 1, True
                        )
                        return
                    else:
                        yield from self.upload_parts(
                            upload, metis_uid, remaining_attempts - 1
                        )
                        return
                raise e

            yield upload

    def create_folder(
        self, project_name: str, bucket_name: str, folder_path: str
    ) -> FoldersAndFilesResponse:
        response = self.session.post(
            self.prepare_url(project_name, "folder", "create", bucket_name, folder_path)
        )
        response_obj = from_json(FoldersAndFilesResponse, response.content)
        return response_obj

    def find(
        self,
        project_name: str,
        bucket_name: str,
        params: List[typing.Mapping],
        limit: Optional[int] = None,
        offset: Optional[int] = None,
        hide_paths=False,
    ) -> FoldersAndFilesResponse:
        limit_set = True

        if not limit:
            limit_set = False
            limit = 1000
        if not offset:
            offset = 0

        args = dict(
            params=params,
            hide_paths=hide_paths,
            limit=limit,
            offset=offset,
        )

        result = FoldersAndFilesResponse(folders=[], files=[])
        while True:
            response = self.session.post(
                self.prepare_url(project_name, "find", bucket_name), json=args
            )
            response_obj = from_json(FoldersAndFilesResponse, response.content)
            result.extend(response_obj)

            if response_obj.empty() or limit_set:
                break

            args["offset"] += max(len(response_obj.files), len(response_obj.folders))

        return result


@serialize
@deserialize
@dataclasses.dataclass
class Attribute:
    attribute_type: str = ""
    link_model_name: Optional[str] = None
    match: Optional[str] = None
    restricted: Optional[bool] = None
    hidden: Optional[bool] = None


@serialize
@deserialize
@dataclasses.dataclass
class Template:
    attributes: Dict[str, Attribute]
    name: str = ""
    identifier: str = ""
    version: int = 0
    parent: str = ""


@serialize
@deserialize
@dataclasses.dataclass
class Model:
    documents: Dict[str, Dict[str, typing.Any]]
    template: Optional[Template] = None
    count: int = 0


@serialize
@deserialize
@dataclasses.dataclass
class RetrievalResponse:
    models: Dict[str, Model]

    def empty(self):
        for model in self.models.values():
            if model.count:
                return False
        return True

    def extend(self, other: "RetrievalResponse"):
        for model_name, model in other.models.items():
            if model_name in self.models:
                self.models[model_name].count += model.count
                self.models[model_name].documents.update(model.documents)
            else:
                self.models[model_name] = model


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

    UPLOAD_PATH_REGEX = re.compile(r"https://.*/([^/]+)/upload/([^/]+)/(.+)")

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

    def read_next_blob(self) -> bytes:
        next_left = self.cur
        next_right = self.cur + self.next_blob_size

        # Advance the position, if we've resumed to a specific place that needs skipping to.
        # We are merely reading data to catch up to our expected left position.
        if next_right < self.read_position:
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


@serialize
@deserialize
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


@serialize
@deserialize
@dataclasses.dataclass
class UpdateRequest:
    revisions: Dict[str, Dict[str, Dict[str, typing.Any]]] = dataclasses.field(
        default_factory=dict
    )
    project_name: str = ""
    dry_run: bool = False

    def extend(self, other: "UpdateRequest"):
        for model_name, docs in other.revisions.items():
            for record_name, revision in docs.items():
                self.update_record(model_name, record_name, revision)

    # TODO: Add more basic validations?
    def validate(self, models: Dict[str, Model]) -> typing.Iterable[str]:
        for model_name, docs in self.revisions.items():
            if model_name not in models:
                yield f"Model '{model_name}' does not exist"
                continue

            template = models[model_name].template
            for record_name, revision in docs.items():
                for attr, value in revision.items():
                    if attr not in template.attributes:
                        yield f"Model '{model_name}', Attribute '{attr}' does not exist"
                        continue

    def update_record(
        self, model_name: str, record_name: str, attrs: Dict[str, typing.Any] = dict()
    ) -> Dict[str, typing.Any]:
        attrs = {k: normalize_for_magma(v) for k, v in attrs.items()}
        record = self.revisions.setdefault(model_name, {}).setdefault(record_name, {})
        record.update(**attrs)
        return record

    def append_table(
        self,
        parent_model_name: str,
        parent_record_name: str,
        model_name: str,
        attrs: Dict[str, typing.Any],
        attribute_name: str,
    ):
        parent_revision = self.update_record(parent_model_name, parent_record_name, {})
        table = parent_revision.setdefault(attribute_name, [])
        id = f"::{model_name}{len(self.revisions.setdefault(model_name, {})) + 1}"
        table.append(id)
        self.update_record(model_name, id, attrs)
        return id


class Magma(EtnaClientBase):
    def retrieve(
        self,
        project_name: str,
        model_name="all",
        attribute_names="all",
        record_names: typing.Union[List[str], "all"] = [],
        page: Optional[int] = None,
        page_size: Optional[int] = None,
        order: Optional[str] = None,
        filter: Optional[str] = None,
        hide_templates=True,
    ) -> RetrievalResponse:
        args = dict(
            project_name=project_name,
            model_name=model_name,
            attribute_names=attribute_names,
            hide_templates=hide_templates,
            record_names=record_names,
        )

        if page is None and record_names:
            page = 1
        if page_size is None and record_names:
            page_size = self.batch_size

        if page is not None:
            args["page"] = page
        if page_size is not None:
            args["page_size"] = page_size

        if order is not None:
            args["order"] = order
        if filter is not None:
            args["filter"] = filter

        response = RetrievalResponse(models={})

        def consume_pages():
            while True:
                try:
                    r = self.session.post(self.prepare_url("retrieve"), json=args)
                except HTTPError as e:
                    if b"not found" in e.response.content:
                        if e.response.status_code == 422:
                            break
                    raise e
                paged = from_json(RetrievalResponse, r.content)
                response.extend(paged)

                if "page" not in args:
                    break

                if paged.empty():
                    break
                args["page"] += 1

        if record_names and isinstance(record_names, list):
            for record_name_batch in batch_iterable(record_names, self.batch_size):
                args["record_names"] = record_name_batch
                args["page"] = 1
                consume_pages()
        else:
            args["record_names"] = record_names
            consume_pages()

        return response

    batch_size = 100

    # TODO: Support batched update, breaking down the request into viable chunks
    def update(self, update: UpdateRequest):
        response = self.session.post(
            self.prepare_url("update"),
            data=to_json(update),
            headers={"Content-Type": "application/json"},
        )
        return from_json(RetrievalResponse, response.content)


def normalize_for_magma(value):
    # Exhaust generators (such as uploads) before consuming their last yielded value
    if isgenerator(value):
        for value in value:
            pass

    if isinstance(value, File):
        return value.as_magma_file_attribute

    if isinstance(value, Upload):
        magma_file_attribute = value.as_magma_file_attribute
        if magma_file_attribute is None:
            raise ValueError(
                f"Could not determine metis path for upload: '{value.upload_path}', possible bug in etna client."
            )
        return magma_file_attribute

    if isinstance(value, Upload):
        return

    if isinstance(value, list):
        return [normalize_for_magma(v) for v in value]

    return value
