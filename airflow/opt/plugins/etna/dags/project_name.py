from airflow import DAG, AirflowException


def project_name_of(dag: DAG):
    parts = dag.owner.split(',')
    if len(parts) == 0:
        raise AirflowException(f"Dag {dag.dag_id} does not have an owner set on any of its tasks, cannot infer project owner.")
    return parts[0]
