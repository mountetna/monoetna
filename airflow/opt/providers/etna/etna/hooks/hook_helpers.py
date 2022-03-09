from airflow.providers.slack.hooks.slack import SlackHook

from etna.dags.project_name import get_project_name
from etna.hooks.connections import find_first_valid_connection
from etna.hooks.etna import EtnaHook
from etna.hooks.git import GitHook


def get_project_slack_hook() -> SlackHook:
    return SlackHook(
        slack_conn_id=find_first_valid_connection(
            f"slack_{get_project_name()}", "slack"
        )
    )


def get_git_hook(remote_path: str) -> GitHook:
    return GitHook(
        git_conn_id=find_first_valid_connection(f"git_{get_project_name()}", "git"),
        remote_path=remote_path,
    )


def get_etna_hook() -> EtnaHook:
    return EtnaHook.for_project()
