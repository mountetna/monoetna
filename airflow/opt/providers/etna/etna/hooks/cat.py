import asyncio
import codecs
import contextlib
from datetime import datetime
import hashlib
import io
import logging
import os
import re
import stat
import subprocess
from paramiko.sftp_client import SFTPClient
from paramiko.sftp_attr import SFTPAttributes
from typing import List, Dict, Optional, ContextManager
from queue import Queue, Empty
from threading  import Thread

import cached_property

from airflow.exceptions import AirflowException
from airflow.models import Connection, Variable
from airflow.operators.python import get_current_context
from airflow.providers.ssh.hooks.ssh import SSHHook
from etna.dags.project_name import get_project_name
from etna.hooks.hook_helpers import RemoteFileBase
from etna.hooks.ssh_base import SSHBase


class CatHook(SSHHook):
    """
    CAT Client Hook to manage the connection
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
            "hidden_fields": ["port", "schema"],
            "relabeling": {},
            "placeholders": {
                "host": "fastq.ucsf.edu",
                "extra": "{\n  \"root_path\": \"SSD\",\n  \"host_key\": \"ssh-<key type> <ssh-key>\"\n}"
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
    def full_path(self) -> str:
        return os.path.join(self.folder_path, self.attrs.filename)

    @property
    def rel_path(self) -> str:
        return self.full_path[1::] if self.full_path[0] == '/' else self.full_path

    def _epoch_s(self, timestamp: datetime):
        return int(timestamp.timestamp())

    def is_in_range(self, batch_start: Optional[datetime] = None, batch_end: Optional[datetime] = None):
        if batch_start is None and batch_end is None:
            return True
        elif batch_start is None:
            return self.mtime <= self._epoch_s(batch_end)
        elif batch_end is None:
            return self.mtime >= self._epoch_s(batch_start)
        else:
            return self._epoch_s(batch_start) <= self.mtime <= self._epoch_s(batch_end)

    @property
    def hash(self):
        # Do NOT use self.atime here, it will update
        #       every time the file is accessed,
        #       which means files always appear newer.
        return f"{self.size}-{self.mtime}"


class Cat(SSHBase):
    variable_root = "cat_ingest_cursor"

    def __init__(self, hook: CatHook):
        self.hook = hook
        self.cursors = {
            "c4": Variable.get(self.variable_key("c4"), default_var={}, deserialize_json=True),
            "metis": Variable.get(self.variable_key("metis"), default_var={}, deserialize_json=True)
        }

    def tail(
            self,
            ignore_directories: List[str],
            magic_string: re.Pattern = re.compile(r".*DSCOLAB.*"),
     ) -> List[SftpEntry]:
        """
        Tails all files that have not been ingested previously.
        If a `magic_string` regex is provided, will only return files that match the provided
            regex.

        params:
          magic_string: regex, the "oligo" we use to identify our projects' files. Default of ".*DSCOLAB.*
          ignore_directories: List[str], list of directories to not scan
        """
        all_files = []

        with self.sftp() as sftp:
            all_files = self._ls_r(sftp, magic_string, ignore_directories)

        return all_files

    def _ls_r(
        self,
        sftp: SFTPClient,
        magic_string: re.Pattern,
        ignore_directories: List[str],
        path: str = None) -> List[SftpEntry]:
        files = []

        if path is None:
            path = self._root_path()

        print(f"Checking path: {path}")

        for entry in sftp.listdir_attr(path):
            sftp_entry = SftpEntry(entry, path)

            if sftp_entry.is_dir() and sftp_entry.name in ignore_directories:
                print(f"{sftp_entry.name} is a skipped dir.")
                continue
            elif sftp_entry.is_dir():
                print(f"{sftp_entry.name} is a directory, going into it")
                files += self._ls_r(sftp, magic_string, ignore_directories, os.path.join(path, sftp_entry.name))
            elif (magic_string.match(sftp_entry.name) and
                 (not self.file_ingested_to_system("c4", sftp_entry) or
                  not self.file_ingested_to_system("metis", sftp_entry))):
                # Re-upload file if the hash has changed (size-modified time-created time).
                print(f"appending {sftp_entry.name} to list of files!")
                files.append(sftp_entry)

        return files

    def _curl_url(self, file: SftpEntry) -> str:
        return f"sftp://{self.hook.connection.host}/{file.full_path}"

    def _curl_authn(self) -> str:
        return f"{self.hook.connection.login}:{self.hook.connection.password}"

    def _curl_hostpubmd5(self) -> str:
        # https://stackoverflow.com/a/66227907
        return hashlib.md5(
            codecs.decode(bytes(self._key_str(),'utf-8'), 'base64')
        ).hexdigest()

    @contextlib.contextmanager
    def retrieve_file(
        self,
        source_file: SftpEntry
    ):
        """
        Opens the given file for download as a file-like object.
        The returned file-like object yields bytes.

        Note:  Ideally, this method is used in combination with 'with' syntax in python, so that the underlying
        data stream is closed after usage.  This is especially performant when code only needs to access a small
        subset of the data.

        eg:
        ```
        with cat.retrieve_file(file) as source_file:
            for line in source_file:
                break
        ```

        params:
          source_file: an SFTP file listing
        """
        def wait_for_proc(pr, qu):
            pr.wait()
            qu.put(pr.returncode)

        rd, wd = os.pipe()

        cmd = [
            "curl",
            "-u",
            self._curl_authn(),
            "--hostpubmd5",
            self._curl_hostpubmd5(),
            "-o",
            "-",
            "-N",
            self._curl_url(source_file)
        ]

        proc = subprocess.Popen(
            cmd,
            stdout=wd,
            stderr=subprocess.PIPE,
            pass_fds=[wd],
            bufsize=0)
        q = Queue()

        closer = Thread(
            target=wait_for_proc,
            args=(proc, q),
            daemon=True
        )

        closer.start()

        os.close(wd)
        try:
            yield open(rd, mode='rb', closefd=False)
            os.close(rd)
        except Exception as e:
            os.close(rd)
            proc.kill()
            raise e

        status = q.get()
        if 0 != status:
            raise AirflowException(f"Failed to run external process, got status code {status}")

    def mark_file_as_ingested(self, ingested_to: str, file: SftpEntry):
        """
        In the cursor, save the fact that the given file's upload was completed.
        """
        if ingested_to not in self.cursors:
            self.cursors[ingested_to] = {}

        self.cursors[ingested_to][file.full_path] = file.hash

    def file_ingested_to_system(self, ingested_to: str, file: SftpEntry) -> bool:
        """
        Returns if the file exists in the system's specific cursor.
        """
        if ingested_to not in self.cursors:
            return False

        return file.full_path in self.cursors[ingested_to] and self.cursors[ingested_to][file.full_path] == file.hash

    def update_cursor(self, postfix: str):
        """
        Save the cursor to the database.
        """
        Variable.set(self.variable_key(postfix), self.cursors[postfix], serialize_json=True)

    def variable_key(self, postfix: str):
        """
        Return the variable key for the current dag.
        """
        try:
            context = get_current_context()

            return f"{self.variable_root}-{postfix}-{context['dag'].dag_id}"
        except AirflowException:
            return f"{self.variable_root}-{postfix}"