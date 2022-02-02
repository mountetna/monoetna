from typing import List

from airflow.exceptions import AirflowNotFoundException, AirflowException
from airflow.hooks.base import BaseHook
from airflow.models import DagRun
from airflow.models.dag import DagStateChangeCallback
from airflow.providers.slack.hooks.slack import SlackHook
from airflow.utils.session import create_session
from airflow.utils.state import State

from airflow.models.taskinstance import Context, TaskInstance

def find_first_valid_connection(*options: str) -> str:
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

def notify_slack_dag_callback(dag_status: str) -> DagStateChangeCallback:
    def state_change_callback(context: Context):
        task_instance: TaskInstance = context["task_instance"]
        conn_id = find_first_valid_connection(
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
