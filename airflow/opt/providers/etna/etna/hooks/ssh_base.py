import base64
import contextlib
import json
import logging
from typing import Tuple
from io import StringIO
import paramiko
from paramiko.sftp_client import SFTPClient

from airflow.exceptions import AirflowException
from airflow.models import Variable
from airflow.operators.python import get_current_context
from airflow.providers.ssh.hooks.ssh import SSHHook

log = logging.getLogger("airflow.task")


class SSHBase(object):
    variable_root = "ssh_ingest_cursor"
    chunk_size: int = 4194304

    def __init__(self, hook: SSHHook):
        self.hook = hook
        self.cursor = Variable.get(self.variable_key, default_var={}, deserialize_json=True)

    @contextlib.contextmanager
    def ssh(self) -> paramiko.SSHClient:
        """
        Configures an SSH connection. Using Python `with` syntax, so that the
        connection is closed after usage.

        eg:
        ```
        with hook.ssh() as ssh:
            ssh.get("/directory")
        ```
        """
        ssh = paramiko.SSHClient()

        keys = ssh.get_host_keys()
        keys.add(
            self.hook.connection.host,
            self._key_type(),
            self._host_key()
        )

        ssh.connect(
            self.hook.connection.host,
            username=self.hook.connection.login,
            password=self.hook.connection.password,
            pkey=self._private_key())

        yield ssh

        ssh.close()

    @contextlib.contextmanager
    def sftp(self) -> SFTPClient:
        """
        Configures an SFTP connection. Using Python `with` syntax, so that the
        connection is closed after usage.

        eg:
        ```
        with hook.sftp() as sftp:
            sftp.get("/directory")
        ```
        """
        with self.ssh() as ssh:
            sftp = ssh.open_sftp()

            yield sftp

            sftp.close()

    def _extra(self) -> dict:
        if (self.hook.connection.extra != ''):
            return json.loads(self.hook.connection.extra)
        return {}

    def _get_extra(self, key, default_value):
        if (key in self._extra() and
            self._extra()[key] != ""):
            return self._extra()[key]

        return default_value

    def _root_path(self) -> str:
        return self._get_extra("root_path", "SSD")

    def _key_components(self) -> Tuple[str, str]:
        host_key_str = self._get_extra('host_key', '')

        if host_key_str == '':
            raise AirflowException("Must provide host_key in the connection's extra options")

        return host_key_str.split(None)[:2]

    def _key_type(self) -> str:
        return self._key_components()[0]

    def _key_str(self) -> str:
        return self._key_components()[1]

    def _host_key(self) -> paramiko.pkey.PKey:
        key_class = {
            "ssh-rsa": paramiko.RSAKey,
            "ssh-ed25519": paramiko.Ed25519Key,
            "ssh-ecdsa": paramiko.ECDSAKey,
            "ssh-dss": paramiko.DSSKey
        }

        if self._key_type() not in key_class:
            raise AirflowException(f"Unsupported SSH key type: {self._key_type()}")

        return key_class[self._key_type()](
            data=base64.b64decode(self._key_str())
        )

    def _private_key(self) -> paramiko.pkey.PKey:
        private_key_str = self._get_extra('private_key', '')

        if '' == private_key_str:
            return None

        return paramiko.rsakey.RSAKey.from_private_key(StringIO(private_key_str))
