import base64
import contextlib
import dataclasses
import json
import typing
from datetime import datetime
from typing import Dict, Optional, List
from urllib.parse import quote

from airflow import DAG
from airflow.operators.python import get_current_context
from serde import serialize, deserialize

import cached_property
import dateutil
import requests
from airflow.exceptions import AirflowException
from airflow.hooks.base import BaseHook
from airflow.models import Connection
from cryptography.hazmat.primitives import serialization, hashes
from cryptography.hazmat.primitives.asymmetric import padding
from requests import HTTPError
from requests.auth import AuthBase
from serde.json import from_json

from etna.dags.project_name import project_name_of
from etna.hooks.keys import prepared_key_from


class EtnaHook(BaseHook):
    """
    Etna Client Hook that can generate task tokens
    """

    conn_name_attr = 'etna_conn_id'
    default_conn_name = 'git_default'
    conn_type = 'git'
    hook_name = 'Git Checkout'
    root_dir = '/opt/airflow/dags/repos'

    remote_path: str
    local_path: str
    etna_conn_id: str
    key: Optional[str]

    # https://github.com/apache/airflow/blob/main/airflow/customized_form_field_behaviours.schema.json
    @staticmethod
    def get_ui_field_behaviour() -> Dict:
        return {
            "hidden_fields": ['port', 'extra'],
            "relabeling": {
                'password': 'RSA key or Janus Token',
                'host': 'environment host postfix',
                'schema': 'rsa or token',
                'login': 'email of janus account (necessary for rsa)'
            },
            'placeholders': {
                'host': '.ucsf.edu',
                'schema': 'token',
                'login': 'you@ucsf.edu'
            }
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

    def get_conn(self) -> Connection:
        return self.get_connection(self.etna_conn_id)

    @cached_property.cached_property
    def connection(self) -> Connection:
        return self.get_conn()


    def get_hostname(self, host_prefix: str):
        return f"{host_prefix}{self.connection.host or '.ucsf.edu'}"

    def get_token_auth(self) -> AuthBase:
        schema = self.connection.schema or 'token'
        if schema == 'rsa':
            with prepared_key_from(self.connection, True) as key_file:
                with self.get_client(None) as session:
                    nonce = Janus(session, self.get_hostname('janus')).nonce()
                return SigAuth(key_file, self.connection.login, nonce)
        elif schema == 'token':
            token = self.connection.password or ''
            return TokenAuth(token.encode('utf8'))
        else:
            raise AirflowException(f"Connection {self.connection.conn_id} has invalid schema: '{self.connection.schema}'")

    def get_task_auth(self, project_name: Optional[str]=None, ready_only=True) -> "TokenAuth":
        if project_name is None:
            dag: DAG = get_current_context()['dag']
            project_name = project_name_of(dag)
        token_auth = self.get_token_auth()
        with self.get_client(token_auth) as session:
            token = Janus(session, self.get_hostname('janus')).generate_token(project_name, ready_only)
        return TokenAuth(token, project_name)

    @contextlib.contextmanager
    def get_client(self, auth: Optional[AuthBase]) -> typing.ContextManager[requests.Session]:
        with requests.Session() as session:
            if auth:
                session.auth = auth
            session.hooks['response'].append(EtnaHook.assert_status)
            yield session

    @contextlib.contextmanager
    def metis(self, project_name: Optional[str] = None, read_only=True) -> typing.ContextManager["Metis"]:
        auth = self.get_task_auth(project_name, read_only)
        with self.get_client(auth) as session:
            yield Metis(session, self.get_hostname('metis'))

    @contextlib.contextmanager
    def janus(self, project_name: Optional[str] = None, read_only=True) -> typing.ContextManager["Janus"]:
        auth = self.get_task_auth(project_name, read_only)
        with self.get_client(auth) as session:
            yield Janus(session, self.get_hostname('janus'))

    @staticmethod
    def assert_status(response: requests.Response, *args, **kwds):
        if 200 <= response.status_code < 300:
            return

        error_message = response.reason

        try:
            err_json = response.json()
            errors = err_json.get('errors', [])
            error = err_json.get('error')

            if error:
                errors.append(error)

            error_message = ', '.join(errors)
            error_message = f"Request failed with {response.status_code}: {error_message}"
        except (TypeError, ValueError, json.JSONDecodeError):
            pass

        raise HTTPError(error_message, response=response)

    def generate_task_token(self):
        pass

class TokenAuth(AuthBase):
    token: bytes
    project_scope: Optional[str]

    def __init__(self, token: bytes, project_scope: Optional[str]=None):
        self.token = token
        self.project_scope = project_scope

    def __call__(self, r: requests.Request) -> requests.Request:
        r.headers['Authorization'] = f'Etna {self.token.decode("ascii")}'
        return r

class SigAuth(AuthBase):
    nonce: bytes
    email: str
    private_key: bytes

    def __init__(self, private_key_file: str, email: str, nonce: bytes):
        self.private_key = open(private_key_file, 'rb').read()
        self.nonce = nonce
        self.email = email

    def __call__(self, r: requests.Request) -> requests.Request:
        key = serialization.load_pem_private_key(self.private_key, None)
        txt_to_sign: bytes = b'.'.join([self.nonce, base64.b64encode(self.email.encode('ascii'))])
        sig = base64.b64encode(key.sign(txt_to_sign, padding.PKCS1v15(
            # If in the future the server changes to PSS
            # mgf=padding.MGF1(hashes.SHA256()),
            # salt_length=padding.PSS.MAX_LENGTH
        ), hashes.SHA256()))
        together: bytes = b'.'.join([txt_to_sign, sig])
        r.headers['Authorization'] = f'Signed-Nonce {together.decode("ascii")}'
        return r

def encode_path(*segments: str) -> str:
    return '/'.join(quote(s) for s in segments)


class EtnaClientBase:
    session: requests.Session
    hostname: str

    def __init__(self, session: requests.Session, hostname: str):
        self.session = session
        self.hostname = hostname

    def prepare_url(self, *path: str):
        return f"https://{self.hostname}/{encode_path(*path)}"

    def get_project_scope(self) -> Optional[str]:
        if isinstance(self.session.auth, TokenAuth):
            return self.session.auth.project_scope
        return None

class Janus(EtnaClientBase):
    def generate_token(self, project_name: str, ready_only=True, token_type='task') -> bytes:
        response = self.session.post(self.prepare_url('api', 'tokens', 'generate'), json=dict(
            token_type=token_type,
            project_name=project_name,
            ready_only=ready_only
        ))
        return response.content

    def nonce(self) -> bytes:
        response = self.session.get(self.prepare_url('api', 'tokens', 'nonce'))
        return response.content

@serialize
@deserialize
@dataclasses.dataclass
class File:
    file_name: str
    size: int
    file_hash: str
    updated_at: str
    file_path: str
    project_name: str
    bucket_name: str
    download_path: str

    @property
    def download_path(self) -> str:
        return encode_path(self.project_name, 'download', self.bucket_name, self.file_path)

    @property
    def updated_at_datetime(self) -> datetime:
        return dateutil.parser.isoparse(self.updated_at)

@serialize
@deserialize
@dataclasses.dataclass
class Folder:
    folder_path: str
    project_name: str
    bucket_name: str
    folder_name: str
    updated_at: str

    @property
    def updated_at_datetime(self) -> datetime:
        return dateutil.parser.isoparse(self.updated_at)

@serialize
@deserialize
@dataclasses.dataclass
class FoldersAndFilesResponse:
    folders: List[Folder]
    files: List[File]

class Metis(EtnaClientBase):
    def list_folder(self, project_name: str, bucket_name: str, folder_path: Optional[str] = None):
        if folder_path:
            response = self.session.get(self.prepare_url(project_name, 'list', bucket_name, folder_path))
        else:
            response = self.session.get(self.prepare_url(project_name, 'list', bucket_name))
        return from_json(FoldersAndFilesResponse, response.content)

    def find(self, project_name: str, bucket_name: str, params: List[typing.Mapping],
             limit: Optional[int] = None, offset: Optional[int] = None, hide_paths=False):
        args = dict(
            params=params,
            hide_paths=hide_paths
        )

        if limit:
            args['limit'] = limit
        if offset is not None:
            args['offset'] = offset

        response = self.session.post(self.prepare_url(project_name, 'find', bucket_name), json=args)
        return from_json(FoldersAndFilesResponse, response.content)



