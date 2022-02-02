from datetime import timedelta, datetime

from airflow import DAG
from airflow.exceptions import AirflowException
from airflow.models import TaskInstance, XCom
from airflow.sensors.base import BaseSensorOperator
from airflow.triggers.temporal import TimeDeltaTrigger
from airflow.utils.session import provide_session
from sqlalchemy.orm import Session

from etna.etls.context import get_batch_range
from etna.xcom.etna_xcom import EtnaDeferredXCom


class AwaitBatches(BaseSensorOperator):
    loader_dag: DAG
    ordering_key: str

    def __init__(self, loader_dag: DAG, ordering_key: str, **kwds):
        super().__init__(**kwds)
        self.loader_dag = loader_dag
        self.ordering_key = ordering_key

    def execute(self, context):
        self.defer(trigger=TimeDeltaTrigger(timedelta(minutes=1)), method_name="check_or_complete")

    @provide_session
    def check_or_complete(self, context, event=None, session: Session=None):
        ti: TaskInstance = context['ti']
        start, end = get_batch_range(context)
        if ti.execution_date < datetime.now() + timedelta(seconds=self.timeout or 60 * 60):
            raise AirflowException(f"Timeout awaiting loaded batch from dag {self.loader_dag.dag_id}")

        # Make sure that we have processed 'past' the current end point, so that our data should be complete.
        row = session.query(XCom).filter(
            XCom.dag_id == self.loader_dag.dag_id,
            XCom.execution_date > end
        ).order_by(XCom.execution_date.asc()).limit(1).first()

        if not row:
            self.defer(trigger=TimeDeltaTrigger(timedelta(minutes=1)), method_name="check_or_complete")

        upper = row.execution_date

        # Make sure that we have processed 'past' the current end point, so that our data should be complete.
        row = session.query(XCom).filter(
            XCom.dag_id == self.loader_dag.dag_id,
            XCom.execution_date <= start
        ).order_by(XCom.execution_date.asc()).limit(1).first()

        if not row:
            self.defer(trigger=TimeDeltaTrigger(timedelta(minutes=1)), method_name="check_or_complete")

        lower = row.execution_date

        return BatchReferenceResult(self.loader_dag.dag_id, self.ordering_key, lower, upper)

class BatchReferenceResult(EtnaDeferredXCom):
    source_dag_id: str
    lower: datetime
    upper: datetime

    def __init__(self, source_dag_id: str, ordering_key: str, lower: datetime, upper: datetime):
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
        return result



