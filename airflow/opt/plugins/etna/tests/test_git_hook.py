from airflow.models import Connection
from unittest import mock

from etna.hooks.git import GitHook
from etna.tests.conftest import NotSoRandom


@mock.patch("tempfile._Random", NotSoRandom)
def test_git_hook_ssh(ssh_git_connection: Connection):
    hook = GitHook(ssh_git_connection.conn_id, 'mountetna/metis.git')
    hook.clone()


@mock.patch("tempfile._Random", NotSoRandom)
def test_git_hook_https(https_git_connection: Connection):
    hook = GitHook(https_git_connection.conn_id, 'mountetna/metis.git')
    hook.clone()
