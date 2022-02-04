import json
from datetime import timedelta, datetime
from typing import Optional, Any

from airflow import DAG
from airflow.exceptions import AirflowException
from airflow.models import TaskInstance, XCom
from airflow.sensors.base import BaseSensorOperator
from airflow.triggers.temporal import TimeDeltaTrigger
from airflow.utils.session import provide_session
from airflow.utils.timezone import utc
from serde.json import from_json
from sqlalchemy.orm import Session

from etna.etls.context import get_batch_range
from etna.xcom.etna_xcom import EtnaDeferredXCom


class AwaitBatches(BaseSensorOperator):
    loader_dag: DAG
    ordering_key: str
    loader_task_id: Optional[str]

    def __init__(self, task_id: str, loader_dag: DAG, element_class: Any, loader_task_id: Optional[str] = None, ordering_key: str = 'updated_at', **kwds):
        super().__init__(task_id=task_id, **kwds)
        self.loader_dag = loader_dag
        self.loader_task_id = loader_task_id
        self.ordering_key = ordering_key
        self.element_class = element_class

    def execute(self, context):
        return self.check_or_complete(context)

    @provide_session
    def check_or_complete(self, context, event=None, session: Session=None):
        ti: TaskInstance = context['ti']
        start, end = get_batch_range(context)

        # Make sure that we have processed 'past' the current end point, so that our data should be complete.
        filters = [
            XCom.dag_id == self.loader_dag.dag_id,
            XCom.execution_date > end
        ]

        if self.loader_task_id:
            filters.append(XCom.task_id == self.loader_task_id)

        row = session.query(XCom).filter(*filters).order_by(XCom.execution_date.asc()).limit(1).first()

        if not row:
            if datetime.utcnow().replace(tzinfo=utc) > (
                    ti.execution_date + timedelta(seconds=self.timeout or 3600)).replace(tzinfo=utc):
                raise AirflowException(f"Timeout awaiting loaded batch from dag {self.loader_dag.dag_id}")
            self.defer(trigger=TimeDeltaTrigger(timedelta(minutes=1)), method_name="check_or_complete")

        upper = row.execution_date

        filters = [
            XCom.dag_id == self.loader_dag.dag_id,
            XCom.execution_date <= start
        ]

        if self.loader_task_id:
            filters.append(XCom.task_id == self.loader_task_id)

        row = session.query(XCom).filter(*filters).order_by(XCom.execution_date.desc()).limit(1).first()

        if not row:
            if datetime.utcnow().replace(tzinfo=utc) > (
                    ti.execution_date + timedelta(seconds=self.timeout or 3600)).replace(tzinfo=utc):
                raise AirflowException(f"Timeout awaiting loaded batch from dag {self.loader_dag.dag_id}")
            self.defer(trigger=TimeDeltaTrigger(timedelta(minutes=1)), method_name="check_or_complete")

        lower = row.execution_date

        return BatchReferenceResult(self.loader_dag.dag_id, self.ordering_key, self.element_class, lower, upper)

class BatchReferenceResult(EtnaDeferredXCom):
    source_dag_id: str
    element_class: Any
    lower: datetime
    upper: datetime

    def __init__(self, source_dag_id: str, ordering_key: str, element_class: Any, lower: datetime, upper: datetime):
        self.source_dag_id = source_dag_id
        self.ordering_key = ordering_key
        self.lower = lower
        self.upper = upper

    @provide_session
    def execute(self, session: Session = None):
        xcoms = session.query(XCom).filter(
            XCom.dag_id == self.source_dag_id,
            XCom.execution_date >= self.lower,
            XCom.execution_date < self.upper,
        ).all()

        result = []
        for xcom in xcoms:
            result.extend(XCom.deserialize_value(xcom))
        result.sort(key=lambda row: row[self.ordering_key])
        return [from_json(self.element_class, json.dumps(row)) for row in result]



