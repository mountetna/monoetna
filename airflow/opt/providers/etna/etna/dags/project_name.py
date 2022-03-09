from typing import Optional

from airflow import DAG, AirflowException
from airflow.operators.python import get_current_context


def get_project_name(dag: Optional[DAG] = None):
    """
    Utility function that retrieves the current project context that a task may be executing in.
    Note that this function may only be called _inside of an executing task_, such as a @link decorated
    function or a @task decorated function.  A DAG's project is determined by the first owner entry, which
    is set by most utilities in the etna library automatically.
    """
    if dag is None:
        dag = get_current_context()["dag"]
    parts = dag.owner.split(",")
    if len(parts) == 0:
        raise AirflowException(
            f"Dag {dag.dag_id} does not have an owner set on any of its tasks, cannot infer project owner."
        )
    return parts[0]
