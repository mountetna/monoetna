"""
`etna` is a collection of airflow related utilities for interacting with
and processing data inside the DSCOLabs Data Library project.

These utilities are meant to supplement those provided directly by the `airflow`
python library, whose basic interfaces and behaviors are compatible and composable
with all the utilities here within.  https://airflow.apache.org/docs/

The top level module exports many common utilities needed by an average etl,
while individual inner packages may contain implementations, interfaces, and
granular functions for nuanced use cases.
"""

from datetime import datetime


def __getattr__(name):
    from .dags.decorators import system_dag, dag, rollup_dag
    from .operators import run_on_docker
    from .etls.decorators import metis_etl
    from .etls.metis import (
        link,
        MetisEtlHelpers,
    )
    from .xcom.etna_xcom import pickled
    from .hooks.hook_helpers import (
        get_project_name,
        get_etna_hook,
        get_git_hook,
        get_project_slack_hook,
    )
    from .hooks.etna import UpdateRequest
    from .operators.rollup_xcom_operator import rollup

    __all__ = list(k for k in locals().keys() if k != "name") + [
        "dags",
        "operators",
        "etls",
        "xcom",
        "hooks",
    ]
    __version__ = "0.0.1"
    __date__ = datetime(2022, 2, 22)

    for k, v in locals().items():
        globals()[k] = v

    if name in locals():
        return locals()[name]
    raise AttributeError(f"{name} does not exist in etna")


# trick the static checker.
if globals().get("notathing", False):
    from .dags.decorators import system_dag, dag, rollup_dag
    from .operators import run_on_docker
    from .etls.decorators import metis_etl
    from .etls.metis import (
        link,
        MetisEtlHelpers,
    )
    from .xcom.etna_xcom import pickled
    from .hooks.hook_helpers import (
        get_project_name,
        get_etna_hook,
        get_git_hook,
        get_project_slack_hook,
    )
    from .hooks.etna import UpdateRequest
    from .operators.rollup_xcom_operator import rollup

    from airflow.models.dag import DAG
    from airflow.models.taskinstance import TaskInstanceKey
    from airflow.utils.task_group import TaskGroup
