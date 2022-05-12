import contextlib
from datetime import datetime, timezone
import io
import logging
import os
import re
from typing import List, Dict, Optional, ContextManager, Tuple, BinaryIO
from ftplib import FTP_TLS, error_reply
from socket import socket

import cached_property
from airflow.hooks.base import BaseHook
from airflow.models import Connection
from etna.dags.project_name import get_project_name


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

    def is_in_range(self, batch_start: Optional[datetime] = None, batch_end: Optional[datetime] = None):
        if batch_start is None and batch_end is None:
            return True
        elif batch_start is None:
            return self.mtime <= batch_end
        elif batch_end is None:
            return self.mtime >= batch_start
        else:
            return batch_start <= self.mtime <= batch_end


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
        """
        Tails all files found in the given `folder_name`, with modified times (`mtime`) in the
        given batch_start and batch_end date range. Recursively searches sub-folders.

        params:
          folder_name: str, the top-level folder to search from.
          batch_start, Optional[datetime], to compare against each file's mtime, via <=
          batch_end, Optional[datetime], to compare against each file's mtime, via >=
        """
        if not Box.valid_folder_name(folder_name):
            raise ValueError(f"Invalid folder name: {folder_name}. Only alphanumeric characters, _, -, and spaces are allowed.")

        with self.ftps() as ftps:
            return self._ls_r(ftps, batch_start=batch_start, batch_end=batch_end)

    def _ls_r(self, ftps: FTP_TLS, path: str = "/", batch_start: Optional[datetime] = None, batch_end: Optional[datetime] = None) -> List[FtpEntry]:
        files = []
        for entry in ftps.mlsd(path):
            ftp_entry = FtpEntry(entry, path)

            if ftp_entry.is_dot():
                continue

            if ftp_entry.is_dir():
                files += self._ls_r(ftps, os.path.join(path, ftp_entry.name))
            elif ftp_entry.is_in_range(batch_start, batch_end):
                files.append(ftp_entry)
        return files

    @contextlib.contextmanager
    def ftps(self) -> FTP_TLS:
        """
        Configures an FTP_TLS connection to Box. Using Pythong `with` syntax, so that the
        connection is closed after usage.

        eg:
        ```
        with box.ftps() as ftps:
            ftps.sendcmd("LIST")
        ```
        """
        ftps = FTP_TLS(self.hook.connection.host)
        ftps.login(user=self.hook.connection.login, passwd=self.hook.connection.password)

        yield ftps

        ftps.quit()

    def retrieve_file(self, ftps: FTP_TLS, file: FtpEntry) -> Tuple[BinaryIO, int]:
        """
        Opens the given file for download into a context as a Tuple(BinaryIO, size).
        The underlying socket object yields bytes objects.

        Note:  Ideally, this method is used in combination with 'with' syntax in python, so that the underlying
        data stream is closed after usage.  This is especially performant when code only needs to access a small
        subset of the data.

        eg:
        ```
        socket = box.retrieve_file(file)
        with socket.makefile('rb') as connection:
            for line in csv.reader(connection):
                break
        ```

        params:
          ftps: an open, FTP_TLS connection
          file: an FTP file listing
        """
        # force binary file transfers via the socket
        ftps.voidcmd('TYPE I')

        # Required because...otherwise the next ftps call may return the status
        #   message from the previous command. :shrug:
        ftps.voidcmd('NOOP')
        return ftps.transfercmd(f"RETR {file.full_path}")

    def remove_file(self, ftps: FTP_TLS, file: FtpEntry):
        """
        Removes the file from the FTP server, if user wants to automatically
        clean up after ingestion to Metis.
        """
        ftps.delete(file.full_path)
