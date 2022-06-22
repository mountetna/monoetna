import contextlib
from datetime import datetime, timezone
import json
import logging
import os
import re
import subprocess
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
                "extra": "{\n  \"root_path\": \"/volume1/SSD\",\n  \"host_key\": \"<ssh-key>\"\n}"
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


class RsyncEntry(RemoteFileBase):
    def __init__(self, raw: str):
        self.raw_line = raw
        self.metadata_delimiter = "+++++++++ "

    @property
    def name(self) -> str:
        return os.path.basename(self.full_path)

    @property
    def size(self) -> int:
        return None

    def is_file(self) -> bool:
        return self.raw_line.startswith(">f")

    def is_dir(self) -> bool:
        return self.raw_line.startswith("cd")

    @property
    def mtime(self) -> datetime:
        return None

    @property
    def full_path(self) -> str:
        return self.raw_line.split(self.metadata_delimiter)[1]

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
        return self.full_path


class Cat(object):
    variable_root = "cat_ingest_cursor"
    ssh_hostfile = "~/.ssh/known_hosts"

    def __init__(self, hook: CatHook):
        self.hook = hook
        self.cursor = Variable.get(self.variable_key, default_var={}, deserialize_json=True)

    def tail(
            self,
            magic_string: re.Pattern = re.compile(r".*DSCOLAB.*"),
     ) -> List[RsyncEntry]:
        """
        Tails all files that have not been ingested previously.
        If a `magic_string` regex is provided, will only return files that match the provided
            regex.

        params:
          magic_string: regex, the "oligo" we use to identify our projects' files. Default of ".*DSCOLAB.*
        """
        self._ensure_cat_host_key()
        log.info("set known_hosts")
        with open(self.ssh_hostfile, 'rb') as host_file:
            for line in host_file.readlines():
                log.info(line)
        all_files = self._rsync()

        return [f for f in all_files if magic_string.match(f.full_path)]

    def _extra(self) -> dict:
        if (self.hook.connection.extra != ''):
            return json.loads(self.hook.connection.extra)
        return {}

    def _root_path(self) -> str:
        if ("root_path" in self._extra() and
            self._extra()["root_path"] != ""):
            return self._extra()["root_path"]

        return "/volume1/SSD"

    def _rsync_connection_str(self) -> str:
        return f"{self.hook.connection.login}@{self.hook.connection.host}:{self._root_path()}"

    def _rsync(self) -> List[RsyncEntry]:
        cmd = [
            "rsync",
            "-avzh",
            "--dry-run",
            "--itemize-changes",
            "--exclude=test",
            "--exclude=Reports",
            "--exclude=Stats",
            self._rsync_connection_str(),
            ".",
        ]

        process = SubprocessHook()

        env: Dict[str, str] = {
            'RSYNC_PASSWORD': self.hook.connection.password
        }

        result = process.run_command(
            cmd,
            env
        )

        if result.exit_code != 0:
            raise AirflowException(
                f"Rsync command failed!  Check logs for more information."
            )

        all_entries = [RsyncEntry(entry) for entry in result.split("\n")]
        return [e for e in all_entries if e.is_file()]

    def _ensure_cat_host_key(self):
        key_exists = self._key_in_hostfile()
        log.info(f"key exists: {key_exists}")
        if not key_exists:
            self._add_key()

    def _key_in_hostfile(self):
        if not os.path.exists(self.ssh_hostfile):
            return False

        with open(self.ssh_hostfile, 'r') as host_file:
            for line in host_file.readlines():
                if (line.startswith(self.hook.connection.host) and
                    -1 != line.find(self._extra()["host_key"])):
                    return True

        return False

    def _add_key(self):
        os.makedirs(os.path.dirname(self.ssh_hostfile), exist_ok=True)
        with open(self.ssh_hostfile, 'a+') as host_file:
            key = self._extra()["host_key"]
            host_file.write(f"{self.hook.connection.host} {key}\n")

        log.info("adding host key")
        log.info(os.path.exists(self.ssh_hostfile))
        with open(self.ssh_hostfile, 'r') as host_file:
            for line in host_file.readlines():
                log.info(line)

    def retrieve_file(self, target_hook: BaseHook, file: RsyncEntry) -> Tuple[BinaryIO, int]:
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

    def mark_file_as_ingested(self, file: RsyncEntry):
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
