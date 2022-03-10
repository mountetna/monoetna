from airflow.plugins_manager import AirflowPlugin

from etna.hooks.git import GitHook


class EtnaPlugin(AirflowPlugin):
    name = "etna_plugin"

    hooks = [GitHook]
