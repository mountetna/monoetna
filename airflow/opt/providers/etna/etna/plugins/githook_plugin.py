from airflow.plugins_manager import AirflowPlugin


from etna.hooks.git import GitHook


class GitHookPlugin(AirflowPlugin):
    name = "githook_plugin"

    hooks = [GitHook]
