from pendulum import datetime
from datetime import timedelta
from typing import Mapping, Union
from typing import Optional

from airflow import DAG
from airflow.models.dag import (
    DagStateChangeCallback,
    ScheduleIntervalArgNotSet,
    ScheduleIntervalArg,
)

from dags.callbacks import notify_slack_dag_callback
from etna.etls.context import batch_start_context_key, batch_end_context_key
from etna.utils.inject import inject

def dag(
        on_failure_callback: Optional[DagStateChangeCallback] = None,
        on_success_callback: Optional[DagStateChangeCallback] = None,
        schedule_interval: ScheduleIntervalArg = ScheduleIntervalArgNotSet,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        inject_params: Mapping[str, str] = {},
        version: Union[int, str] = "",
        **kwds,
):
    """
    Creates a dag by running the wrapped function in a dag context whose dag_id will be the function's name,
    and whose description will be its doc string.
    :return: Dah'ya like dags?
    """

    def instantiate_dag(fn):
        with DAG(
                dag_id=fn.__name__ + str(version),
                description=fn.__doc__,
                on_failure_callback=on_failure_callback,
                on_success_callback=on_success_callback,
                start_date=start_date,
                end_date=end_date,
                schedule_interval=schedule_interval,
                **kwds,
        ) as dag:
            inject(fn, inject_params)
        return dag

    return instantiate_dag


def etl(project_name: str, start_date: datetime, interval: timedelta, version: Union[int, str], inject_params: Mapping = {}):
    def instantiate_dag(fn):
        return dag(
            start_date=start_date,
            schedule_interval=interval,
            catchup=True,
            inject_params=dict(
                batch_start_date="{{ " + batch_start_context_key + " }}",
                batch_end_date="{{ " + batch_end_context_key + " }}",
                version=version,
                **inject_params
            ),
            default_args=dict(
                owner=project_name,
                retries=5,
            ),
            version=version,
        )(fn)

    return instantiate_dag

system_epoch = datetime(2021, 12, 22, 16, 56, 3, 185905)

# A dag context
def system_dag(interval: timedelta):
    def instantiate_dag(fn):
        return dag(
            start_date=system_epoch,
            schedule_interval=interval,
            default_args=dict(
                owner='administration',
                retries=5,
            ),
            # default_args=dict(
            on_failure_callback=notify_slack_dag_callback("failed: "),
            # ),
            catchup=False,
        )(fn)

    return instantiate_dag
