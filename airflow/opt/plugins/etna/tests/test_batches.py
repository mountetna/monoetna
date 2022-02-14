from datetime import datetime

from airflow import DAG
from airflow.models import TaskInstance
from airflow.operators.dummy import DummyOperator

from etna.etls.batches import maybe_set_task_lower_bound, get_task_lower_bound, LOWEST_BOUND


def test_maybe_set_task_lower_bound():
    dag = DAG('dag_id')
    task = DummyOperator(task_id='task_id', dag=dag, start_date=datetime(2020, 1, 1))
    context = dict(ti=TaskInstance(task=task))

    assert get_task_lower_bound(context) == LOWEST_BOUND

    maybe_set_task_lower_bound(context, datetime(2020, 1, 1))
    maybe_set_task_lower_bound(context, datetime(2019, 1, 1))
    assert get_task_lower_bound(context) == datetime(2020, 1, 1)
