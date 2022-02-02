import base64
import contextlib
import json
from typing import Dict, Optional

import os

import cached_property
import requests
import rsa
from airflow.hooks.base import BaseHook
from airflow.models import Connection
from requests import HTTPError
from requests.auth import AuthBase


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
            "hidden_fields": ['port', 'extra', 'schema', 'login'],
            "relabeling": {
                'password': 'long lived token to generate task tokens from',
                'host': 'host part after service name (eg .ucsf.edu, -stage.ucsf.edu, .development.local)'
            },
            'placeholders': {
                'host': '.ucsf.edu',
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

    def __init__(self, etna_conn_id: str, remote_path: str, local_path: Optional[str] = None) -> None:
        super().__init__()
        self.etna_conn_id = etna_conn_id
        self.remote_path = remote_path

        if local_path is None:
            local_path = os.path.basename(self.remote_path)

        if not local_path.startswith('/'):
            local_path = os.path.join(self.root_dir, local_path)

        self.local_path = local_path
        self.key = None

    def get_conn(self) -> Connection:
        return self.get_connection(self.etna_conn_id)

    @cached_property.cached_property
    def connection(self) -> Connection:
        return self.get_conn()

    @contextlib.contextmanager
    def get_client(self, token: str) -> requests.Session:
        with requests.Session() as session:
            session.auth = TokenAuth(token)
            session.hooks['response'].append(EtnaHook.assert_status)
            yield session

    @contextlib.contextmanager
    def janus(self, token: str) -> "Janus":
        with self.get_client(token) as session:
            yield Janus(session)

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
    def __init__(self, token):
        self.token = token

    def __call__(self, r: requests.Request) -> requests.Request:
        r.headers['Authorization'] = f'Etna {self.token}'
        return r

class NonceAuth(AuthBase):
    def __init__(self, private_key_file: str, email: str, nonce: bytes):
        self.private_key_file = private_key_file
        self.nonce = nonce
        self.email = email

    def __call__(self, r: requests.Request) -> requests.Request:
        pk = rsa.PrivateKey.load_pkcs1(open(self.private_key_file, 'rb').read())
        txt_to_sign: bytes = b'.'.join([self.nonce, base64.b64encode(self.email.encode('utf8'))])
        sig = base64.b64encode(rsa.sign(txt_to_sign, pk, ' SHA-256'))
        together: bytes = b'.'.join([txt_to_sign, sig])
        r.headers['Authorization'] = f'Signed-Nonce {together.decode("ascii")}'
        return r

class Janus:
    session: requests.Session

    def __init__(self, session: requests.Session):
        self.session = session

    def generate_token(self, ready_only = True, token_type = 'task'):
        pass

