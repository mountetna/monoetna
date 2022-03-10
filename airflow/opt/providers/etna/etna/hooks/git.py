import contextlib
import tempfile
from typing import Dict, Optional

import os

import stat

import cached_property
from airflow import AirflowException
from airflow.hooks.base import BaseHook
from airflow.hooks.subprocess import SubprocessHook
from airflow.models import Connection

from etna.hooks.keys import prepared_key_from


class GitHook(BaseHook):
    """
    Hook to connect to a remote git repository using an ssh key.
    """

    conn_name_attr = "git_conn_id"
    default_conn_name = "git_default"
    conn_type = "git"
    hook_name = "Git Checkout"
    root_dir = "/opt/airflow/dags/repos"

    remote_path: str
    local_path: str
    git_conn_id: str
    key: Optional[str]

    # https://github.com/apache/airflow/blob/main/airflow/customized_form_field_behaviours.schema.json
    @staticmethod
    def get_ui_field_behaviour() -> Dict:
        return {
            "hidden_fields": ["port", "extra"],
            "relabeling": {
                "schema": "ssh or https",
                "password": "http password, or ssh private key contents",
                "login": "http user, or ssh user (git for github)",
            },
            "placeholders": {
                "host": "github.com",
                "schema": "ssh",
                "login": "git",
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

    def __init__(
        self, git_conn_id: str, remote_path: str, local_path: Optional[str] = None
    ) -> None:
        super().__init__()
        self.git_conn_id = git_conn_id
        self.remote_path = remote_path

        if local_path is None:
            local_path = os.path.basename(self.remote_path)

        if not local_path.startswith("/"):
            local_path = os.path.join(self.root_dir, local_path)

        self.local_path = local_path
        self.key = None

    @contextlib.contextmanager
    def key_path(self) -> Optional[str]:
        if self.connection.schema != "https":
            with prepared_key_from(self.connection) as file_name:
                yield file_name
        else:
            yield None

    def get_conn(self) -> Connection:
        return self.get_connection(self.git_conn_id)

    @cached_property.cached_property
    def connection(self) -> Connection:
        return self.get_conn()

    def run_git(self, working_dir: str, *command: str):
        with self.key_path() as key_file_path:
            process = SubprocessHook()
            env: Dict[str, str] = dict(**os.environ)

            with tempfile.NamedTemporaryFile(delete=True) as ask_pass_file:
                if key_file_path is not None:
                    env["GIT_SSH_COMMAND"] = f"ssh -i {key_file_path} -F /dev/null"
                elif self.connection.password:
                    os.chmod(ask_pass_file.name, stat.S_IRUSR | stat.S_IXUSR)
                    ask_pass_file.write("!#/usr/bin/env bash")
                    ask_pass_file.write(f"exec echo {self.connection.password}")
                    env["GIT_ASKPASS"] = ask_pass_file.name

                result = process.run_command(
                    ["git"] + list(command), env, cwd=working_dir
                )

            if result.exit_code != 0:
                raise AirflowException(
                    f"git command failed!  Check logs for more information."
                )

    def clone(self):
        os.makedirs(os.path.dirname(self.local_path), exist_ok=True)

        if self.connection.schema == "https":
            host = self.connection.host or "github.com"
            user = self.connection.login
            if user:
                clone_url = f"https://{user}@{host}/{self.remote_path}"
            else:
                clone_url = f"https://{host}/{self.remote_path}"
        else:
            login = self.connection.login or "git"
            host = self.connection.host or "github.com"
            clone_url = f"{login}@{host}:{self.remote_path}"

        return self.run_git("/opt/airflow", "clone", clone_url, self.local_path)

    def pull(self, *args: str):
        return self.run_git(self.local_path, "pull", *args)

    def checkout(self, *args: str):
        return self.run_git(self.local_path, "checkout", *args)

    def commit(self, *args: str):
        return self.run_git(self.local_path, "commit", *args)

    def push(self, *args: str):
        return self.run_git(self.local_path, "push", *args)
