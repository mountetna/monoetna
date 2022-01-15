import base64
import io
import json
import tarfile
from datetime import datetime, timedelta
from os.path import basename
from typing import Callable
from typing import Optional, List, Any

from airflow import DAG
from airflow.exceptions import AirflowNotFoundException, AirflowException
from airflow.hooks.base import BaseHook
from airflow.models import DagRun
from airflow.models.dag import (
    DagStateChangeCallback,
    ScheduleIntervalArgNotSet,
    ScheduleIntervalArg,
)
from airflow.models.taskinstance import Context, TaskInstance
from airflow.providers.slack.hooks.slack import SlackHook
from airflow.utils.session import create_session
from airflow.utils.state import State

from etna.operators.docker_operator import DockerOperator
from etna.operators.docker_operator_base import SwarmSharedData
from etna.operators.swarm_operator import DockerSwarmOperator
from etna.get_provider_info import get_provider_info

system_epoch = datetime(2021, 12, 22, 16, 56, 3, 185905)


def _find_first_valid_connection(*options: str) -> str:
    last_e = None
    for option in options:
        try:
            BaseHook.get_connection(option)
            return option
        except AirflowNotFoundException as e:
            last_e = e
    if last_e:
        raise last_e

    raise AirflowNotFoundException("No connection option present!")


def _notify_slack_dag_callback(dag_status: str) -> DagStateChangeCallback:
    def state_change_callback(context: Context):
        task_instance: TaskInstance = context["task_instance"]
        conn_id = _find_first_valid_connection(
            f"slack_notifications_{task_instance.dag_id}", "slack_notifications"
        )

        with create_session() as session:
            task_instance = session.merge(task_instance)
            dag_run: DagRun = task_instance.dag_run
            task_instances: List[TaskInstance] = list(dag_run.task_instances)
            failed_instances: List[TaskInstance] = [
                ti for ti in task_instances if ti.state in State.failed_states
            ]

        reason = context.get("reason", "")

        failed_message = ""
        if failed_instances:
            failed_message = f"\nFailed task logs: {'  '.join(ti.log_url for ti in failed_instances)}"

        slack = SlackHook(slack_conn_id=conn_id)
        conn = slack.get_connection(conn_id)

        if not getattr(conn, "host", None):
            raise AirflowException("Missing channel(host) in Slack connection")

        api_params = {
            "channel": conn.host,
            "username": conn.login or "Airflow",
            "text": f"Dag {task_instance.dag_id} {dag_status}{reason}!{failed_message}",
            "icon_url": conn.schema
            or "https://raw.githubusercontent.com/apache/airflow/main/airflow/www/static/pin_100.png",
        }

        slack.call("chat.postMessage", json=api_params)

    return state_change_callback


def system_dag(interval: timedelta):
    def instantiate_dag(fn):
        return dag(
            start_date=system_epoch,
            schedule_interval=interval,
            # default_args=dict(
            on_failure_callback=_notify_slack_dag_callback("failed: "),
            # ),
            catchup=False,
        )(fn)

    return instantiate_dag


def dag(
    on_failure_callback: Optional[DagStateChangeCallback] = None,
    on_success_callback: Optional[DagStateChangeCallback] = None,
    schedule_interval: ScheduleIntervalArg = ScheduleIntervalArgNotSet,
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
    **kwds,
):
    """
    Creates a dag by running the wrapped function in a dag context whose dag_id will be the function's name,
    and whose description will be its doc string.
    :return: Dah'ya like dags?
    """

    def instantiate_dag(fn):
        with DAG(
            dag_id=fn.__name__,
            description=fn.__doc__,
            on_failure_callback=on_failure_callback,
            on_success_callback=on_success_callback,
            start_date=start_date,
            end_date=end_date,
            schedule_interval=schedule_interval,
            **kwds,
        ) as dag:
            fn()
        return dag

    return instantiate_dag


def run_in_container(
    task_id: str,
    source_container: str,
    command: List[str],
    input_bytes: Optional[bytes] = None,
    use_compose: bool = True,
    input_file: Optional[str] = None,
    input_dir: Optional[str] = None,
    output_json: bool = False,
    output_b64: bool = False,
):
    # Uses the containers created by the local development docker-compose system.
    if use_compose:
        source_container = f"monoetna_{source_container}_1"
    serialize_last_output, swarm_shared_data = _prepare_input_outputs(
        input_bytes, input_dir, input_file, output_b64, output_json
    )

    return DockerOperator(
        task_id=task_id,
        source_container_name=source_container,
        command=command,
        swarm_shared_data=swarm_shared_data,
        serialize_last_output=serialize_last_output,
    )


def run_in_swarm(
    task_id: str,
    source_service: str,
    command: List[str],
    input_bytes: Optional[bytes] = None,
    input_file: Optional[str] = None,
    input_dir: Optional[str] = None,
    include_external_network: bool = False,
    output_json: bool = False,
    output_b64: bool = False,
):
    serialize_last_output, swarm_shared_data = _prepare_input_outputs(
        input_bytes, input_dir, input_file, output_b64, output_json
    )

    return DockerSwarmOperator(
        task_id=task_id,
        source_service=source_service,
        command=command,
        include_external_networks=include_external_network,
        swarm_shared_data=swarm_shared_data,
        serialize_last_output=serialize_last_output,
    )


def _prepare_input_outputs(input_bytes, input_dir, input_file, output_b64, output_json):
    swarm_shared_data: List[SwarmSharedData] = []
    serialize_last_output: Optional[Callable[[bytes], Any]] = None
    if input_bytes:
        swarm_shared_data.append(
            SwarmSharedData(
                data=input_bytes,
                remote_path="/input",
            )
        )
    elif input_file:
        swarm_shared_data.append(
            SwarmSharedData(
                data=open(input_file, "r").read(),
                remote_path="/input/" + basename(input_file),
            )
        )
    elif input_dir:
        buff = io.BytesIO(b"")
        tar = tarfile.open(fileobj=buff)
        tar.add(input_dir, arcname=basename(input_dir))
        tar.close()

        swarm_shared_data.append(
            SwarmSharedData(
                data=buff.read(),
                remote_path="/input/" + basename(input_dir) + ".tar",
            )
        )
    if output_json:
        serialize_last_output = json.loads
    elif output_b64:
        serialize_last_output = base64.b64decode
    return serialize_last_output, swarm_shared_data
