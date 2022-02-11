from datetime import timedelta, datetime
from typing import Optional, Any, Tuple, Generic, TypeVar, List

import dateutil.parser
from airflow.exceptions import AirflowRescheduleException, AirflowException
from airflow.models import TaskInstance, XCom, BaseOperator, DagRun, Variable, XCOM_RETURN_KEY
from airflow.models.taskinstance import Context
from airflow.models.xcom_arg import XComArg
from airflow.utils.session import provide_session
from airflow.utils.timezone import utc
from sqlalchemy.orm import Session

from etna.xcom.etna_xcom import EtnaDeferredXCom

T = TypeVar("T")

batch_start_context_key = 'prev_data_interval_end_success'
batch_end_context_key = 'data_interval_end'

# A theoretical floor for date ranges.
LOWEST_BOUND = datetime(2010, 1, 1).replace(tzinfo=utc)

# Used to determine where this particular task's batch processing was last able to successfully complete.
# Tasks added into the system after a dag has been run will use this value to help backfill themselves on dagruns
# they did not participate in.
def _get_task_lower_bound(context: Context):
    ti: TaskInstance = context['ti']
    lower = Variable.get(f"{ti.dag_id}-{ti.task_id}-lower-bound", LOWEST_BOUND.isoformat())
    return dateutil.parser.isoparse(lower)

# After, and only after, a task completes within a given batch range, should this be given.
def _set_task_lower_bound(context: Context, end: datetime):
    ti: TaskInstance = context['ti']
    Variable.set(f"{ti.dag_id}-{ti.task_id}-lower-bound", end.isoformat())

def _get_dag_batch_range(context: Context) -> Tuple[datetime, datetime]:
    ti: TaskInstance = context['ti']
    start_date = context[batch_start_context_key] or LOWEST_BOUND
    end_date = context[batch_end_context_key]
    return start_date, end_date

def get_batch_range(context: Context) -> Tuple[datetime, datetime]:
    return _get_dag_batch_range(context)

@provide_session
def _get_batch_range(context: Context, session: Session = None) -> Tuple[datetime, datetime]:
    """
    Selects a processing time range such that if this task is 'behind' the dag run, the time range will instead
    include a past time range based on the last known successful processing according to _get_task_lower_bound
    (A Variable) and the dag run that follows it.  Net result is that the date range will assist in backfilling a task.

    One should only use this function in conjunction with BackfillingTaskOperator that should wrap another, existing
    operator.  In general, most etl related decorators and operators do this by default, check their documentation
    to ensure.

    Using this function in an operator that does not correctly invoke _set_task_lower_bound would result in
    processing a dead range over and over.
    """
    # Tasks can be added to a dag after the fact, but need to perform catchup broader than the dag's own lower bound.
    # We keep a running variable of task specific lower bounds in order to manage this.
    task_lower_bound: datetime = _get_task_lower_bound(context)
    dag_lower, dag_upper = _get_dag_batch_range(context)

    # If our task state is not behind the dagrun, use the dag run's configuration.
    if task_lower_bound >= dag_lower:
        return dag_lower, dag_upper

    upper_dr = session.query(DagRun).filter(DagRun.execution_date > task_lower_bound).order_by(DagRun.execution_date.asc()).first()
    if upper_dr is None:
        task_upper_bound = dag_upper
        # dagrun time boundaries should overlap slightly to ensure precision isn't a factor in missed entries.
        task_upper_bound += timedelta(seconds=1)
    else:
        task_upper_bound = upper_dr.execution_date

    task_upper_bound = task_upper_bound.replace(microsecond=0)
    task_lower_bound = task_lower_bound.replace(microsecond=0)

    return task_lower_bound, task_upper_bound

def enable_task_backfill(op: BaseOperator):
    """
    Generally not needed to be called directly, used by the etl decorator.

    Mutates a given operator by replacing its execute and xcom functionality to support automated task backfill.
    execute will reschedule itself to collect and store xcom results over all batch ranges,
    xcom_pull and xcom_push will support selecting over these batch ranges when necessary.
    """
    if isinstance(op, XComArg):
        op = XComArg.operator

    def xcom_pull(
        context: Any,
        task_ids: Optional[List[str]] = None,
        dag_id: Optional[str] = None,
        key: str = XCOM_RETURN_KEY,
        include_prior_dates: Optional[bool] = None,
    ):
        ti: TaskInstance = context['ti']

        if dag_id is None:
            dag_id = ti.dag_id

        if not task_ids:
            task_ids = [ti.task_id]

        start, stop = _get_batch_range(context)

        return [BatchReferenceResult(dag_id, task_id, start, stop).execute() for task_id in task_ids]

    def xcom_push(
            context: Any,
            key: str,
            value: Any,
            execution_date: Optional[datetime] = None,
    ):
        ti: TaskInstance = context['ti']

        if execution_date is None:
            execution_date, _ = _get_batch_range(context)

        XCom.set(
            key=key,
            value=value,
            task_id=ti.task_id,
            dag_id=ti.dag_id,
            execution_date=execution_date,
        )

    def execute(context: Context):
        ti: TaskInstance = context['ti']
        start, stop = _get_batch_range(context)
        dag_start, dag_stop = _get_dag_batch_range(context)

        try:
            context[batch_start_context_key] = start
            context[batch_end_context_key] = stop

            result = ti.task.__class__.execute(ti.task, context)
        finally:
            context[batch_start_context_key] = dag_start
            context[batch_end_context_key] = dag_stop

        if op.do_xcom_push:
            xcom_push(context, XCOM_RETURN_KEY, result, start)
        _set_task_lower_bound(context, stop)

        if stop < dag_stop:
            # Keep reprocessing to catchup
            raise AirflowRescheduleException(datetime.utcnow().replace(tzinfo=utc))

        return result

    op.xcom_push = xcom_push
    op.xcom_pull = xcom_pull
    op.execute = execute

class BatchReferenceResult(EtnaDeferredXCom, Generic[T]):
    source_dag_id: str
    source_task_id: str
    lower: datetime
    upper: datetime

    batch_ordering_key = 'updated_at'

    def __init__(self, source_dag_id: str, source_task_id: str, lower: datetime, upper: datetime):
        self.source_dag_id = source_dag_id
        self.source_task_id = source_task_id
        self.lower = lower
        self.upper = upper

    def _find_xcom_bound(self, session, filter, asc) -> Optional[datetime]:
        order = XCom.execution_date.asc()
        if not asc:
            order = XCom.execution_date.desc()

        xcom = session.query(XCom).filter(
            XCom.dag_id == self.source_dag_id,
            XCom.task_id == self.source_task_id,
            filter
        ).order_by(order).first()

        if xcom:
            return xcom.execution_date

    @provide_session
    def execute(self, session: Session = None) -> List[T]:
        """
        Attempts to find the batched dataset belonging to to XComs sourced from the dask_id and task_id given
        in the constructor that belongs within the time range provided by upper and lower (incl both ends).
        To this end, it tries to select XCom entries that belong within the 'best' bounds that can be determined
        from available data.  That data is then joined together and returned.
        """
        lower = self._find_xcom_bound(session, XCom.execution_date <= self.lower, False) or \
                self._find_xcom_bound(session, XCom.execution_date > self.lower, True)
        upper = self._find_xcom_bound(session, XCom.execution_date >= self.upper, True) or \
                self._find_xcom_bound(session, XCom.execution_date < self.upper, False)

        if lower is None or upper is None:
            return []

        xcoms = session.query(XCom).filter(
            XCom.dag_id == self.source_dag_id,
            XCom.task_id == self.source_task_id,
            XCom.execution_date >= lower,
            XCom.execution_date <= upper,
        ).with_entities(XCom.value).order_by(XCom.execution_date).all()

        result = []
        for xcom in xcoms:
            next_in = XCom.deserialize_value(xcom)
            if not isinstance(next_in, list):
                raise AirflowException(f"tasks inside etl should return None or a list, found {type(next_in)} from {self.source_task_id}/{self.source_task_id}")
            result.extend(next_in)
        return result
