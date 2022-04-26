# https://airflow.apache.org/docs/apache-airflow-providers/
from etna.metrics.rollup_metrics import register_rollup_metrics


def get_provider_info():
    register_rollup_metrics()

    return {
        "package-name": "etna",
        "name": "etna",
        "description": "",
        "connection-types": [
            {
                "connection-type": "git",
                "hook-class-name": "etna.hooks.git.GitHook",
            },
            {
                "connection-type": "etna",
                "hook-class-name": "etna.hooks.etna.EtnaHook",
            },
            {
                "connection-type": "box",
                "hook-class-name": "etna.hooks.box.BoxHook",
            },
        ],
    }
