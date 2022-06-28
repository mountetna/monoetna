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

from mountetna import Janus, Magma, Metis, TokenAuth, SigAuth

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
                nonce = Janus(None, self.get_hostname("janus")).nonce()
                return SigAuth(key_file, self.connection.login, nonce)
        elif schema == "token":
            token = self.connection.password or ""
            return TokenAuth(token)
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
        # Unfortunately, 'read only' here becomes a viewer token, which does not have
        # read permission to many buckets.  Unfortunately, read_only does not mean 'can read',
        # it means 'cannot write, but also likely cannot read.'
        # We are forced then to always create read_only=False tokens for now.
        token = Janus(token_auth, self.get_hostname("janus")).generate_token(
            project_name, False
        )
        return TokenAuth(token, project_name)

    @contextlib.contextmanager
    def metis(
        self, project_name: Optional[str] = None, read_only=True
    ) -> typing.ContextManager["Metis"]:
        auth = self.get_task_auth(project_name, read_only)
        yield Metis(auth, self.get_hostname("metis"))

    @contextlib.contextmanager
    def janus(
        self, project_name: Optional[str] = None, read_only=True
    ) -> typing.ContextManager["Janus"]:
        auth = self.get_task_auth(project_name, read_only)
        yield Janus(auth, self.get_hostname("janus"))

    @contextlib.contextmanager
    def magma(
        self, project_name: Optional[str] = None, read_only=True
    ) -> typing.ContextManager["Magma"]:
        auth = self.get_task_auth(project_name, read_only)
        yield Magma(auth, self.get_hostname("magma"))


    def generate_task_token(self):
        pass
