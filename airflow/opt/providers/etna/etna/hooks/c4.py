import contextlib
import logging
import os
import re
from typing import Dict, Optional, ContextManager

import cached_property
from paramiko.sftp_client import SFTPClient

from airflow.exceptions import AirflowException
from airflow.models import Connection
from airflow.providers.ssh.hooks.ssh import SSHHook
from etna.dags.project_name import get_project_name
from etna.hooks.ssh_base import SSHBase
from etna.hooks.cat import CatHook, SftpEntry


class C4Hook(SSHHook):
    """
    C4 Client Hook to manage the connection
    """

    conn_name_attr = "c4_conn_id"
    default_conn_name = "c4_default"
    conn_type = "c4"
    hook_name = "C4 Connection"
    root_dir = "/opt/airflow/dags/repos"

    remote_path: str
    local_path: str
    c4_conn_id: str
    key: Optional[str]

    # https://github.com/apache/airflow/blob/main/airflow/customized_form_field_behaviours.schema.json
    @staticmethod
    def get_ui_field_behaviour() -> Dict:
        return {
            "hidden_fields": ["port", "schema"],
            "relabeling": {},
            "placeholders": {
                "host": "c4.ucsf.edu",
                "extra": "{\n  \"root_path\": \"/airflow/ingest\",\n  \"host_key\": \"ssh-<key type> <ssh-key>\"\n}",
                "private_key": "private-key-string"
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

    def __init__(self, c4_conn_id: str) -> None:
        super().__init__(
            ssh_conn_id=c4_conn_id,
            remote_host=self.get_connection(c4_conn_id).host
        )
        self.c4_conn_id = c4_conn_id

    @classmethod
    def for_project(cls, project_name: Optional[str] = None):
        if project_name is None:
            project_name = get_project_name()
        return cls(f"c4_{project_name}")

    def get_conn(self) -> Connection:
        return self.get_connection(self.c4_conn_id)

    @cached_property.cached_property
    def connection(self) -> Connection:
        return self.get_conn()

    @contextlib.contextmanager
    def c4(
        self
    ) -> ContextManager["C4"]:
        yield C4(self)


log = logging.getLogger("airflow.task")


class C4(SSHBase):
    chunk_size: int = 4194304

    def __init__(self, hook: C4Hook):
        self.hook = hook

    def _alphanumeric(self, param: str, include_period: bool = False) -> bool:
        if include_period:
            safe_regex = re.compile(r"^[\w\.]+$")
        else:
            safe_regex = re.compile(r"^\w+$")
        return bool(safe_regex.match(param))

    def _validate_cat_config(self, cat_hook: CatHook):
        if (not self._alphanumeric(cat_hook.connection.login) or
            not self._alphanumeric(cat_hook.connection.password) or
            not self._alphanumeric(cat_hook.connection.host, True)):
            raise AirflowException("Invalid CAT connection configuration. Please make sure all fields are alphanumeric.")

    def upload_file(
        self,
        sftp: SFTPClient,
        folder_path: str,
        file_name: str,
        source_file: SftpEntry,
        cat_hook: CatHook
     ):
        """
        Stream the file to C4, to the given path.

        params:
          sftp: a connected SFTPClient instance
          folder_path: str, the destination path on C4 (relative to the hook's configured root_path)
          file_name: str, the file name to save to, on C4
          source_file: SftpEntry, the source file to download
          cat_hook: CatHook, the CAT hook with SFTP login information
        """
        # Make sure the cat config is safe
        self._validate_cat_config(cat_hook)

        # Ensure the folder path exists on C4
        full_path = os.path.join(self._root_path(), folder_path)
        traversed = []
        for folder in full_path.split("/"):
            if folder == "":
                continue

            path_to_folder = os.path.join("/", "/".join(traversed), folder)
            try:
                sftp.chdir(path_to_folder)
            except IOError:
                sftp.mkdir(path_to_folder)
            traversed.append(folder)

        cmd = f"lftp sftp://{cat_hook.connection.login}:{cat_hook.connection.password}@{cat_hook.connection.host}  -e \"set xfer:clobber on; get {source_file.full_path} -o {os.path.join(full_path, file_name)}; bye\""

        with self.ssh() as ssh, ssh.get_transport() as transport, transport.open_channel("session") as channel:
            channel.set_combine_stderr(True)

            channel.exec_command(cmd)

            exit_status = channel.recv_exit_status()
            if channel.recv_ready():
                data = channel.recv(1024)
                while data:
                    print(data)
                    data = channel.recv(1024)

            if 0 == exit_status:
                print("LFTP GET completed.")
            else:
                raise AirflowException(f"Error reading from channel")
