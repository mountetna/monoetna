import json
import math
from bisect import bisect_left, bisect_right
from datetime import timedelta, datetime
import dateutil.parser
from typing import Optional, Any, Tuple

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
    batch_ordering_key: str
    loader_task_id: Optional[str]

    def __init__(self, task_id: str, loader_dag: DAG, loader_task_id: Optional[str] = None, batch_ordering_key: str = 'updated_at', **kwds):
        super().__init__(task_id=task_id, **kwds)
        self.loader_dag = loader_dag
        self.loader_task_id = loader_task_id
        self.batch_ordering_key = batch_ordering_key

    def execute(self, context):
        return self.check_or_complete(context)

    def _select_xcom(self, session: Session, cond, asc: bool):
        filters = [
            XCom.dag_id == self.loader_dag.dag_id,
            cond
        ]

        if self.loader_task_id:
            filters.append(XCom.task_id == self.loader_task_id)

        ord_by = XCom.execution_date
        if asc:
            ord_by = ord_by.asc()
        else:
            ord_by = ord_by.desc()

        return session.query(XCom).filter(*filters).order_by(ord_by).limit(1).first()

    def _defer_or_timeout(self, ti: TaskInstance):
        if datetime.utcnow().replace(tzinfo=utc) > (
                ti.execution_date + timedelta(seconds=self.timeout or 3600)).replace(tzinfo=utc):
            raise AirflowException(f"Timeout awaiting loaded batch from dag {self.loader_dag.dag_id}")
        self.defer(trigger=TimeDeltaTrigger(timedelta(minutes=1)), method_name="check_or_complete")

    @provide_session
    def check_or_complete(self, context, event=None, session: Session=None):
        ti: TaskInstance = context['ti']
        start, end = get_batch_range(context)
        abs = session.query(XCom).all()

        # Make sure that we have processed 'past' the current end point, so that our data should be complete.
        row = self._select_xcom(session, XCom.execution_date > end, asc=True)
        if not row:
            self._defer_or_timeout(ti)

        upper = row.execution_date

        # Find the bound closest to the lower end.  It usually is that we have a loading batch that ran
        # before us, but in initial load cases that may not be.
        row = self._select_xcom(session, XCom.execution_date <= start, asc=False)
        if not row:
            row = self._select_xcom(session, XCom.execution_date >= start, asc=True)
        if not row:
            self._defer_or_timeout(ti)

        lower = row.execution_date

        return BatchReferenceResult(self.loader_dag.dag_id, self.batch_ordering_key, (lower, start), (upper, end))

class BatchReferenceResult(EtnaDeferredXCom):
    source_dag_id: str
    lower: Tuple[datetime, datetime]
    upper: Tuple[datetime, datetime]

    def __init__(self, source_dag_id: str, batch_ordering_key: str, lower: Tuple[datetime, datetime], upper: Tuple[datetime, datetime]):
        self.source_dag_id = source_dag_id
        self.batch_ordering_key = batch_ordering_key
        self.lower = lower
        self.upper = upper

    @provide_session
    def execute(self, session: Session = None):
        execution_lower, batch_lower = self.lower
        execution_upper, batch_upper = self.upper

        xcoms = session.query(XCom).filter(
            XCom.dag_id == self.source_dag_id,
            XCom.execution_date >= execution_lower,
            XCom.execution_date < execution_upper,
        ).all()

        result = []
        for xcom in xcoms:
            result.extend([(dateutil.parser.isoparse(r[self.batch_ordering_key]), id(r), r) for r in XCom.deserialize_value(xcom)])
        result.sort(key=lambda row: row[1])
        left = bisect_left(result, (batch_lower, 0, {}))
        right = bisect_right(result, (batch_upper, math.inf, {}))
        return [r[2] for r in result[left:right]]



