import base64
import contextlib
from datetime import datetime, timezone
import json
import logging
import os
import re
import paramiko
import stat
from paramiko.sftp_client import SFTPClient
from paramiko.sftp_attr import SFTPAttributes
from typing import List, Dict, Optional, ContextManager, Tuple, BinaryIO

import cached_property

from airflow.exceptions import AirflowException
from airflow.hooks.base import BaseHook
from airflow.hooks.subprocess import SubprocessHook
from airflow.models import Connection, Variable
from airflow.operators.python import get_current_context
from airflow.providers.ssh.hooks.ssh import SSHHook
from etna.dags.project_name import get_project_name
from etna.hooks.hook_helpers import RemoteFileBase


class CatHook(SSHHook):
    """
    Cat Client Hook to manage the connection
    """

    conn_name_attr = "cat_conn_id"
    default_conn_name = "cat_default"
    conn_type = "cat"
    hook_name = "CAT Connection"
    root_dir = "/opt/airflow/dags/repos"

    remote_path: str
    local_path: str
    cat_conn_id: str
    key: Optional[str]

    # https://github.com/apache/airflow/blob/main/airflow/customized_form_field_behaviours.schema.json
    @staticmethod
    def get_ui_field_behaviour() -> Dict:
        return {
            "hidden_fields": ["port"],
            "relabeling": {},
            "placeholders": {
                "host": "fastq.ucsf.edu",
                "schema": "ssh",
                "extra": "{\n  \"root_path\": \"/volume1/SSD\",\n  \"host_key\": \"<ssh-key>\",\n  \"key_type\": \"<ssh-rsa, ssh-dss, etc.\"\n}"
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

    def __init__(self, cat_conn_id: str) -> None:
        super().__init__(
            ssh_conn_id=cat_conn_id,
            remote_host=self.get_connection(cat_conn_id).host
        )
        self.cat_conn_id = cat_conn_id

    @classmethod
    def for_project(cls, project_name: Optional[str] = None):
        if project_name is None:
            project_name = get_project_name()
        return cls(f"cat_{project_name}")

    def get_conn(self) -> Connection:
        return self.get_connection(self.cat_conn_id)

    @cached_property.cached_property
    def connection(self) -> Connection:
        return self.get_conn()

    @contextlib.contextmanager
    def cat(
        self
    ) -> ContextManager["Cat"]:
        yield Cat(self)


log = logging.getLogger("airflow.task")


class SftpEntry(RemoteFileBase):
    def __init__(self, attrs: SFTPAttributes, folder_path: str):
        self.attrs = attrs
        self.folder_path = folder_path

    @property
    def name(self) -> str:
        return self.attrs.filename

    @property
    def size(self) -> int:
        return self.attrs.st_size

    def is_file(self) -> bool:
        return stat.S_ISREG(self.attrs.st_mode)

    def is_dir(self) -> bool:
        return stat.S_ISDIR(self.attrs.st_mode)

    @property
    def mtime(self) -> datetime:
        return self.attrs.st_mtime

    @property
    def atime(self) -> datetime:
        return self.attrs.st_atime

    @property
    def full_path(self) -> str:
        return os.path.join(self.folder_path, self.attrs.filename)

    @property
    def rel_path(self) -> str:
        return self.full_path[1::] if self.full_path[0] == '/' else self.full_path

    def is_in_range(self, batch_start: Optional[datetime] = None, batch_end: Optional[datetime] = None):
        if batch_start is None and batch_end is None:
            return True
        elif batch_start is None:
            return self.mtime <= batch_end
        elif batch_end is None:
            return self.mtime >= batch_start
        else:
            return batch_start <= self.mtime <= batch_end

    @property
    def hash(self):
        return f"{self.size}-{self.mtime}-{self.atime}"


class Cat(object):
    variable_root = "cat_ingest_cursor"
    ssh_hostfile = "~/.ssh/known_hosts"

    def __init__(self, hook: CatHook):
        self.hook = hook
        self.cursor = Variable.get(self.variable_key, default_var={}, deserialize_json=True)

    def tail(
            self,
            magic_string: re.Pattern = re.compile(r".*DSCOLAB.*"),
     ) -> List[SftpEntry]:
        """
        Tails all files that have not been ingested previously.
        If a `magic_string` regex is provided, will only return files that match the provided
            regex.

        params:
          magic_string: regex, the "oligo" we use to identify our projects' files. Default of ".*DSCOLAB.*
        """
        self._ensure_ssh_known_hosts()

        all_files = []

        with self.sftp() as sftp:
            all_files = self._ls_r(sftp)

        return [f for f in all_files if magic_string.match(f.full_path)]

    @contextlib.contextmanager
    def sftp(self) -> SFTPClient:
        """
        Configures an SFTP connection to Box. Using Python `with` syntax, so that the
        connection is closed after usage.

        eg:
        ```
        with cat.sftp() as sftp:
            sftp.get("/directory")
        ```
        """
        ssh = paramiko.SSHClient()

        keys = ssh.get_host_keys()
        keys.add(
            self.hook.connection.host,
            self._key_type(),
            self._host_key()
        )
        keys.save(self.ssh_hostfile)

        ssh.connect(
            self.hook.connection.host,
            username=self.hook.connection.login,
            password=self.hook.connection.password)

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

    def _key_type(self) -> str:
        return self._get_extra("key_path", "ssh-rsa")

    def _host_key(self) -> str:
        return self._get_extra("host_key", "")

    def _ls_r(self, sftp: SFTPClient, path: str = None) -> List[SftpEntry]:
        files = []

        if path is None:
            path = self._root_path()

        for entry in sftp.listdir_iter(path):
            print(entry)
            sftp_entry = SftpEntry(entry, path)

            if sftp_entry.is_dir():
                files += self._ls_r(sftp, os.path.join(path, sftp_entry.name))
            elif sftp_entry.full_path not in self.cursor or self.cursor[sftp_entry.full_path] != sftp_entry.hash:
                # Re-upload file if the hash has changed (size-modified time-created time).
                files.append(sftp_entry)

        return files

    def _ensure_ssh_known_hosts(self):
        if not os.path.exists(self.ssh_hostfile):
            os.makedirs(os.path.dirname(self.ssh_hostfile), exist_ok=True)

    def retrieve_file(self, target_hook: BaseHook, file: SftpEntry) -> Tuple[BinaryIO, int]:
        """
        Opens the given file for download into a context as a Tuple(BinaryIO, size).
        The underlying socket object yields bytes objects.

        Note:  Ideally, this method is used in combination with 'with' syntax in python, so that the underlying
        data stream is closed after usage.  This is especially performant when code only needs to access a small
        subset of the data.

        eg:
        ```
        socket = cat.retrieve_file(file)
        with socket.makefile('rb') as connection:
            for line in csv.reader(connection):
                break
        ```

        params:
          target_hook: Hook with a connection to the target system, where to save the file
          file: an Rsync file listing
        """
        # Will have to use different command per target system
        #       for C4, is lget
        #       for Metis, just stream to it?
        pass

    def mark_file_as_ingested(self, file: SftpEntry):
        """
        In the cursor, save the fact that the given file's upload was completed.
        """
        self.cursor[file.full_path] = file.hash

    def update_cursor(self):
        """
        Save the cursor to the database.
        """
        Variable.set(self.variable_key, self.cursor, serialize_json=True)

    @property
    def variable_key(self):
        """
        Return the variable key for the current dag.
        """
        try:
            context = get_current_context()

            return f"{self.variable_root}-{context['dag'].dag_id}"
        except AirflowException:
            return self.variable_root
