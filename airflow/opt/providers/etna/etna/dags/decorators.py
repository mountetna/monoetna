from pendulum import datetime
from datetime import timedelta
from typing import Mapping, Union, List, Callable
from typing import Optional

from airflow import DAG
from airflow.models.dag import (
    DagStateChangeCallback,
    ScheduleIntervalArgNotSet,
    ScheduleIntervalArg,
)

from etna.dags.callbacks import notify_slack_dag_callback
from etna.utils.inject import inject

system_epoch = datetime(2000, 1, 1)


def dag(
    # A callback function that is invoked when a dag fails.  This can be used to send notifications.
    on_failure_callback: Optional[DagStateChangeCallback] = None,
    # A callback function that is invoked when a dag succeeds.  This can be used to send notifications.
    on_success_callback: Optional[DagStateChangeCallback] = None,
    # Generally a timedelta value, or a cron like string specification, indicating the desired interval of execution
    schedule_interval: ScheduleIntervalArg = ScheduleIntervalArgNotSet,
    # The point of the desired first execution.  Generally, this should be set to some recent date that marks the start of a project.
    start_date: Optional[datetime] = None,
    # Optional 'end point' of processing, after which no further runs will be scheduled, however manual runs can always be triggered.
    end_date: Optional[datetime] = None,
    # A mapping of values to parameter names that will be injected into the decorated function.
    inject_params: Mapping[str, str] = {},
    # A value that is appended to the name of the resulting DAG, allowing one to easily force reprocessing of a DAG when
    # significant changes to its implementation occur.
    version: Union[int, str] = "",
    **kwds,
):
    """
    Simple convenience decorator that wraps airflow.dag decorator with some basic defaults and the ability
    to 'inject' parameters into the callable.  Generally not used directly, use other more specific decorators
    such as system_dag or metis_etl.
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


def system_dag(interval: timedelta) -> Callable[[Callable], DAG]:
    """
    Decorator that converts the decorator function into a DAG by the same name
    using any tasks instantiated within the function.  For use with 'system'
    level dags whose owner should be 'administration'.
    """

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
