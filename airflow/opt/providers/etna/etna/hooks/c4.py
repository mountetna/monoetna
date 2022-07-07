import base64
import contextlib
from datetime import datetime, timezone
import io
import json
import logging
import os
import re
import time
import paramiko
import stat
from paramiko.sftp_client import SFTPClient
from paramiko.sftp_attr import SFTPAttributes
from typing import List, Dict, Optional, ContextManager, Tuple, BinaryIO

import cached_property

from airflow.exceptions import AirflowException
from airflow.models import Connection, Variable
from airflow.operators.python import get_current_context
from airflow.providers.ssh.hooks.ssh import SSHHook
from etna.dags.project_name import get_project_name
from etna.hooks.hook_helpers import RemoteFileBase
from etna.hooks.ssh_base import SSHBase
from etna.hooks.cat import SftpEntry


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
                "extra": "{\n  \"root_path\": \"/airflow/ingest\",\n  \"host_key\": \"ssh-<key type> <ssh-key>\"\n}"
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
    variable_root = "c4_ingest_cursor"
    chunk_size: int = 4194304

    def __init__(self, hook: C4Hook):
        self.hook = hook
        self.cursor = Variable.get(self.variable_key, default_var={}, deserialize_json=True)

    def upload_file(
        self,
        folder_path: str,
        file_name: str,
        file_obj: io.BytesIO,
        file_size: int
     ):
        """
        Stream the file to C4, to the given path.

        params:
          folder_path: str, the destination path on C4 (relative to the hook's configured root_path)
          file_name: str, the file name to save to, on C4
          file_obj: io, the file-like object that is a data source
          file_size: int, size of the file in bytes
        """
        with self.sftp() as sftp:
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

            with sftp.file(os.path.join(full_path, file_name), "w") as target_file:
                bytes_read = 0
                should_log = True
                while bytes_read < file_size:
                    chunk = file_obj.read(self.chunk_size)
                    target_file.write(chunk)
                    bytes_read += len(chunk)
                    # Only log every 5 seconds, to save log space...
                    time_check = int(time.time())
                    if time_check % 5 == 0 and should_log:
                        log.info("Upload progress: " + "{:.2%}".format(bytes_read / file_size))
                        should_log = False
                    elif time_check % 5 != 0 and not should_log:
                        should_log = True
