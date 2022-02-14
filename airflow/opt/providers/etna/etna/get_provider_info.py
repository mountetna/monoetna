# https://airflow.apache.org/docs/apache-airflow-providers/
def get_provider_info():
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
            }
        ]
    }