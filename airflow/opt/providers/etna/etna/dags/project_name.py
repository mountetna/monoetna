from typing import Optional

from airflow import DAG, AirflowException
from airflow.operators.python import get_current_context


def get_project_name(dag: Optional[DAG] = None):
    if dag is None:
        dag = get_current_context()["dag"]
    parts = dag.owner.split(",")
    if len(parts) == 0:
        raise AirflowException(
            f"Dag {dag.dag_id} does not have an owner set on any of its tasks, cannot infer project owner."
        )
    return parts[0]
