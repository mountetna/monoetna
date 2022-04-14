import contextlib
import dataclasses
from datetime import datetime
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


@serialize
@deserialize
@dataclasses.dataclass
class BoxFile:
    file_name: str = ""
    folder_path: Optional[int] = None

    @classmethod
    def is_file(cls, line):
        return line.startswith("get -O")

    @classmethod
    def from_lftp_line(cls, host, folder_name, line):
        sub_path = line.split(os.path.join(host, "Root", folder_name, ""))[1]
        path, filename = os.path.split(sub_path)
        return from_json(BoxFile, json.dumps({
            "file_name": filename,
            "folder_path": os.path.join(folder_name, path)
        }))

    def __str__(self):
        return self.full_path

    @property
    def path(self):
        return os.path.join(self.folder_path, self.file_name)

    @property
    def full_path(self):
        return os.path.join("/", "Root", self.folder_path, self.file_name)


@serialize
@deserialize
@dataclasses.dataclass
class BoxFilesResponse:
    files: List[BoxFile] = dataclasses.field(default_factory=list)

    def empty(self):
        return not self.files

    def extend(self, other: "BoxFilesResponse"):
        self.files.extend(other.files)



log = logging.getLogger("airflow.task")

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
     ) -> List[BoxFile]:
        if not Box.valid_folder_name(folder_name):
            raise ValueError(f"Invalid folder name: {folder_name}. Only alphanumeric characters, _, -, and spaces are allowed.")

        if batch_start and batch_end:
            args = dict(
                batch_start=batch_start.date().isoformat(),
                batch_end=batch_end.date().isoformat(),
                folder_name=folder_name,
            )
        else:
            args = dict(
                folder_name=folder_name,
            )

        response = self._lftp(**args)

        results: List[BoxFile] = []

        for line in response.split("\n"):
            if line and BoxFile.is_file(line):
                results.append(BoxFile.from_lftp_line(self.hook.connection.host, folder_name, line))

        return results

    def _lftp(self, folder_name: str, batch_start: Optional[datetime] = None, batch_end: Optional[datetime] = None):
        """
        Constructs the LFTP command to get a listing of files with a timestamp in the given date range.
        """
        command = [
            "lftp",
            "-u",
            f"{self.hook.connection.login},{self.hook.connection.password}",
            f"{self.hook.connection.schema}://{self.hook.connection.host}",
            "-e",
            f"{self._mirror_cmd(folder_name=folder_name, batch_start=batch_start, batch_end=batch_end)}",
        ]

        with contextlib.ExitStack() as stack:
            tmpdir = stack.enter_context(tempfile.TemporaryDirectory(prefix='box_lftp'))

            process = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=tmpdir,
                encoding='utf-8'
            )

            if process.returncode != 0:
                raise AirflowException(
                    f"LFTP command failed!\n{process.stderr}"
                )

        return process.stdout

    def _mirror_cmd(self, folder_name: str, batch_start: Optional[datetime] = None, batch_end: Optional[datetime] = None):
        command = [
            "mirror",
            "--dry-run",
        ]

        if batch_start:
            command.append(f"--newer-than={batch_start}")

        if batch_end:
            command.append(f"--older-than={batch_end}")

        command.extend([
            f"/Root/{folder_name}/",
            ".;",
            "bye"
        ])

        return ' '.join(command)

    @contextlib.contextmanager
    def ftps(self) -> FTP_TLS:
        """
        Configures an FTP_TLS connection to Box.
        """
        ftps = FTP_TLS(self.hook.connection.host)
        ftps.login(user=self.hook.connection.login, passwd=self.hook.connection.password)

        yield ftps

        ftps.quit()

    def file_size(self, ftps: FTP_TLS, file: BoxFile) -> int:
        """
        Returns the size of the given file.
        """
        ftps.cwd(os.path.dirname(file.full_path))
        return ftps.size(file.file_name)

    def retrieve_file(self, ftps: FTP_TLS, file: BoxFile) -> io.BufferedReader:
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
        ftps.cwd(os.path.dirname(file.full_path))

        io_obj = io.BytesIO()
        ftps.retrbinary(f"RETR {file.file_name}", io_obj.write)
        io_obj.seek(0)
        return io_obj
