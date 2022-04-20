import contextlib
import dataclasses
from datetime import datetime, timezone
import dateutil
import io
import json
import logging
import os
import re
import subprocess
import tempfile
from typing import List, Dict, Optional, ContextManager
from ftplib import FTP_TLS

import cached_property
from airflow import AirflowException
from airflow.hooks.base import BaseHook
from airflow.models import Connection
from etna.dags.project_name import get_project_name
from serde import serialize, deserialize
from serde.json import from_json, to_json

from etna.utils.streaming import iterable_to_stream


class BoxHook(BaseHook):
    """
    Box Client Hook to manage the connection
    """

    conn_name_attr = "box_conn_id"
    default_conn_name = "box_default"
    conn_type = "box"
    hook_name = "Box Connection"
    root_dir = "/opt/airflow/dags/repos"

    remote_path: str
    local_path: str
    box_conn_id: str
    key: Optional[str]

    # https://github.com/apache/airflow/blob/main/airflow/customized_form_field_behaviours.schema.json
    @staticmethod
    def get_ui_field_behaviour() -> Dict:
        return {
            "hidden_fields": ["port", "extra"],
            "relabeling": {},
            "placeholders": {
                "host": "ftp.box.com",
                "schema": "ftps"
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

    def __init__(self, box_conn_id: str) -> None:
        super().__init__()
        self.box_conn_id = box_conn_id

    @classmethod
    def for_project(cls, project_name: Optional[str] = None):
        if project_name is None:
            project_name = get_project_name()
        return cls(f"box_{project_name}")

    def get_conn(self) -> Connection:
        return self.get_connection(self.box_conn_id)

    @cached_property.cached_property
    def connection(self) -> Connection:
        return self.get_conn()

    @contextlib.contextmanager
    def box(
        self
    ) -> ContextManager["Box"]:
        yield Box(self)


log = logging.getLogger("airflow.task")


class FtpEntry(object):
    def __init__(self, tuple, parent_path: str):
        self.name = tuple[0]
        self.metadata = tuple[1]
        self.folder_path = parent_path

    @property
    def size(self) -> int:
        return int(self.metadata['size'])

    def is_file(self) -> bool:
        return self.metadata['type'] == 'file'

    def is_dir(self) -> bool:
        return self.metadata['type'] == 'dir'

    @property
    def mtime(self) -> datetime:
        return datetime.strptime(self.metadata['modify'], "%Y%m%d%H%M%S.000").replace(tzinfo=timezone.utc)

    def is_dot(self) -> bool:
        return self.name == "." or self.name == ".."

    @property
    def full_path(self) -> str:
        return os.path.join(self.folder_path, self.name)

    @property
    def rel_path(self) -> str:
        return self.full_path[1::] if self.full_path[0] == '/' else self.full_path

class Box(object):
    def __init__(self, hook: BoxHook):
        self.hook = hook

    @classmethod
    def valid_folder_name(cls, folder_name: str):
        regex = re.compile(r"^[\w -]+$")
        return regex.match(folder_name)

    def tail(
            self,
            folder_name: str,
            batch_start: Optional[datetime] = None,
            batch_end: Optional[datetime] = None,
     ) -> List[FtpEntry]:
        if not Box.valid_folder_name(folder_name):
            raise ValueError(f"Invalid folder name: {folder_name}. Only alphanumeric characters, _, -, and spaces are allowed.")

        with self.ftps() as ftps:
            all_files = self._ls_r(ftps)

            return [f for f in all_files if self._is_in_range(f, batch_start, batch_end)]

    def _ls_r(self, ftps: FTP_TLS, path: str = "/") -> List[FtpEntry]:
        files = []
        for entry in ftps.mlsd(path):
            ftp_entry = FtpEntry(entry, path)

            if ftp_entry.is_dot():
                continue

            if ftp_entry.is_dir():
                files += self._ls_r(ftps, os.path.join(path, ftp_entry.name))
            else:
                files.append(ftp_entry)
        return files

    def _is_in_range(self, file: FtpEntry, batch_start: Optional[datetime] = None, batch_end: Optional[datetime] = None):
        if batch_start is None and batch_end is None:
            return True
        elif batch_start is None:
            return file.mtime <= batch_end
        elif batch_end is None:
            return file.mtime >= batch_start
        else:
            return batch_start <= file.mtime <= batch_end

    @contextlib.contextmanager
    def ftps(self) -> FTP_TLS:
        """
        Configures an FTP_TLS connection to Box.
        """
        ftps = FTP_TLS(self.hook.connection.host)
        ftps.login(user=self.hook.connection.login, passwd=self.hook.connection.password)

        yield ftps

        ftps.quit()

    def retrieve_file(self, ftps: FTP_TLS, file: FtpEntry) -> io.BufferedReader:
        """
        Opens the given file for download into a context as a python io (file like) object.
        The underlying io object yields bytes objects.

        Note:  Ideally, this method is used in combination with 'with' syntax in python, so that the underlying
        data stream is closed after usage.  This is especially performant when code only needs to access a small
        subset of the data.

        eg:
        ```
        with box.open_file(file) as open_file:
          for line in csv.reader(open_file):
             break
        ```
        """
        ftps.cwd(file.folder_path)

        io_obj = io.BytesIO()
        ftps.retrbinary(f"RETR {file.name}", io_obj.write)
        io_obj.seek(0)
        return io_obj

    def remove_file(self, ftps: FTP_TLS, file: FtpEntry):
        """
        Removes the file from the FTP server. This is so we can keep track of which files
        have been ingested, since we can't rely on upload / file modified timestamps.
        """
        ftps.cwd(file.folder_path)

        ftps.delete(file.name)
