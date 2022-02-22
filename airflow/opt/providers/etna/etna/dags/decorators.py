from pendulum import datetime
from datetime import timedelta
from typing import Mapping, Union, List
from typing import Optional

from airflow import DAG
from airflow.models.dag import (
    DagStateChangeCallback,
    ScheduleIntervalArgNotSet,
    ScheduleIntervalArg,
)

from etna.dags.callbacks import notify_slack_dag_callback
from etna.utils.inject import inject

system_epoch = datetime(2021, 12, 22, 16, 56, 3, 185905)


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


def system_dag(interval: timedelta):
    def instantiate_dag(fn):
        return dag(
            start_date=system_epoch,
            schedule_interval=interval,
            default_args=dict(
                owner="administration",
                retries=3,
            ),
            # default_args=dict(
            # on_failure_callback=notify_slack_dag_callback("failed: "),
            # ),
            catchup=False,
        )(fn)

    return instantiate_dag
