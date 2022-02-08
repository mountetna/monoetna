import json
import math
from bisect import bisect_left, bisect_right
from datetime import timedelta, datetime
import dateutil.parser
from typing import Optional, Any, Tuple, Generic, TypeVar, List, Union

from airflow import DAG
from airflow.decorators import task
from airflow.exceptions import AirflowException, AirflowRescheduleException
from airflow.models import TaskInstance, XCom, BaseOperator, DagRun
from airflow.models.taskinstance import Context
from airflow.models.xcom_arg import XComArg
from airflow.operators.python import get_current_context
from airflow.sensors.base import BaseSensorOperator
from airflow.triggers.temporal import TimeDeltaTrigger
from airflow.utils.session import provide_session
from airflow.utils.timezone import utc
from serde.json import from_json
from sqlalchemy.orm import Session

from etna.etls.context import get_batch_range
from etna.xcom.etna_xcom import EtnaDeferredXCom, pickled

T = TypeVar("T")

@provide_session
def expand_full_batch(
        source_dag_id: str,
        source_task_id: str,
        loader_batch: List[T], # Type inference help
        session: Optional[Session]=None,
) -> Optional["BatchReferenceResult[T]"]:
    context = get_current_context()
    batch_start, batch_end = get_batch_range(context)
    synchronous_with_dag = source_dag_id == context['ti'].dag_id

    def _select_xcom(cond, asc: bool):
        # In this case, we can guaranteed that the schedule of these sources align, and we do not need to defer
        # until a discrete future execution to ensure full batch consumption.

        filters = [
            XCom.dag_id == source_dag_id,
            XCom.task_id == source_task_id,
            cond
        ]

        ord_by = XCom.execution_date
        if asc:
            ord_by = ord_by.asc()
        else:
            ord_by = ord_by.desc()

        return session.query(XCom).filter(*filters).order_by(ord_by).limit(1).first()

    if synchronous_with_dag:
        row = _select_xcom(XCom.execution_date >= batch_end, asc=True)
    else:
        row = _select_xcom(XCom.execution_date > batch_end, asc=True)

    if not row:
        return None
    upper = row.execution_date

    # Find the bound closest to the lower end.  It usually is that we have a loading batch that ran
    # before us, but in initial load cases that may not be.
    row = _select_xcom(XCom.execution_date <= batch_start, asc=False)
    if not row:
        row = _select_xcom(XCom.execution_date >= batch_start, asc=True)
    if not row:
        return None

    lower = row.execution_date

    return BatchReferenceResult(source_dag_id, 'updated_at', (lower, batch_start), (upper, batch_end))

def task(operator: BaseOperator):
    return task

class BackfillingTaskOperator(BaseOperator):
    def __init__(self, orig: BaseOperator):
        self.__dict__['_orig'] = orig

    def __getattr__(self, item):
        return getattr(self._orig, item)

    def __setattr__(self, item, value):
        return setattr(self._orig, item, value)

    @provide_session
    def execute(self, context: Context, session: Session = None):
        dag: DAG = context['dag']
        ti: TaskInstance = context['ti']
        # session.query(DagRun).filter(DagRun.execution_date < ti.execution_date).order_by(DagRun.execution_date.asc())

        raise AirflowRescheduleException('')


class AwaitBatches(BaseSensorOperator):
    loader_dag: DAG
    batch_ordering_key: str
    loader_task_id: str

    def __init__(self, task_id: str, loader_dag_or_id: Union[DAG, str], loader_task_or_id: Union[BaseOperator, str],
                 batch_ordering_key: str = 'updated_at', **kwds):
        super().__init__(task_id=task_id, **kwds)
        self.loader_dag = loader_dag_or_id.dag_id if isinstance(loader_dag_or_id, DAG) else loader_dag_or_id
        self.loader_task_id = loader_task_or_id.task_id if isinstance(loader_task_or_id, BaseOperator) else loader_task_or_id
        self.batch_ordering_key = batch_ordering_key

    def execute(self, context):
        return self.check_or_complete(context)

    def _defer_or_timeout(self, ti: TaskInstance):
        if datetime.utcnow().replace(tzinfo=utc) > (
                ti.execution_date + timedelta(seconds=self.timeout or 3600)).replace(tzinfo=utc):
            raise AirflowException(f"Timeout awaiting loaded batch from dag {self.loader_dag.dag_id}")
        self.defer(trigger=TimeDeltaTrigger(timedelta(minutes=1)), method_name="check_or_complete")

    @provide_session
    def check_or_complete(self, context, event=None, session: Session = None) -> "BatchReferenceResult[Any]":
        result = expand_full_batch(self.loader_dag.dag_id, self.loader_task_id)
        if result is None:
            self._defer_or_timeout(context['ti'])
        return result


class BatchReferenceResult(EtnaDeferredXCom, Generic[T]):
    source_dag_id: str
    lower: Tuple[datetime, datetime]
    upper: Tuple[datetime, datetime]

    def __init__(self, source_dag_id: str, batch_ordering_key: str, lower: Tuple[datetime, datetime],
                 upper: Tuple[datetime, datetime]):
        self.source_dag_id = source_dag_id
        self.batch_ordering_key = batch_ordering_key
        self.lower = lower
        self.upper = upper

    @provide_session
    def execute(self, session: Session = None) -> List[T]:
        execution_lower, batch_lower = self.lower
        execution_upper, batch_upper = self.upper

        xcoms = session.query(XCom).filter(
            XCom.dag_id == self.source_dag_id,
            XCom.execution_date >= execution_lower,
            XCom.execution_date < execution_upper,
        ).all()

        result = []
        for xcom in xcoms:
            result.extend([(dateutil.parser.isoparse(
                r[self.batch_ordering_key] if isinstance(r, dict) else getattr(r, self.batch_ordering_key)), id(r), r)
                           for r in XCom.deserialize_value(xcom)])
        result.sort(key=lambda row: row[1])
        left = bisect_left(result, (batch_lower, 0, {}))
        right = bisect_right(result, (batch_upper, math.inf, {}))
        return pickled([r[2] for r in result[left:right]])
