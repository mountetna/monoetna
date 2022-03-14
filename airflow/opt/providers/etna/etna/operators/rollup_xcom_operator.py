from datetime import datetime
from functools import reduce
from typing import Any, Optional, Callable, Iterable, Tuple, Dict

from airflow.models import BaseOperator, XCom, TaskInstance
from airflow.models.xcom_arg import XComArg
from airflow.utils.session import provide_session
from airflow.utils.state import State

from etna.etls.etl_task_batching import get_batch_range
from sqlalchemy.orm import Session

from etna.metrics.rollup_metrics import Rollup


class RollupState(Rollup):
    last_execution_date: datetime
    result: Rollup

    def __init__(self, last_execution_date: datetime, result: Rollup):
        self.result = result
        self.last_execution_date = last_execution_date

    def measure(self) -> Iterable[Tuple[Dict[str, str], int]]:
        return self.result.measure()

class rollup:
    concat: Callable[[Rollup, Rollup], Rollup]

    def __init__(self, concat: Callable[[Rollup, Rollup], Rollup]):
        self.concat = concat

    def __call__(self, arg: XComArg):
        op = RollupXcomOperator(concat=self.concat, arg=arg)
        arg >> op
        return op

class RollupXcomOperator(BaseOperator):
    concat: Callable[[Any, Any], Any]
    target_dag_id: str
    target_task_id: str

    def __init__(self,
                 concat: Callable[[Rollup, Rollup], Rollup],
                 arg: Optional[XComArg] = None,
                 target_dag_id: Optional[str] = None,
                 target_task_id: Optional[str] = None,
                 task_id: Optional[str] = None, **kwargs
                 ):
        if arg is not None:
            target_dag_id = target_dag_id or arg.operator.dag_id
            target_task_id = target_task_id or arg.operator.task_id

        super().__init__(task_id=task_id or f"rollup_{target_dag_id}_{target_task_id}", **kwargs)
        self.concat = concat

        self.target_dag_id = target_dag_id
        self.target_task_id = target_task_id

    @provide_session
    def execute(self, context: Any, session: Session = None):
        start, stop = get_batch_range(context)

        result: Optional[RollupState] = None

        ti: TaskInstance = context['ti']
        previous: Optional[TaskInstance] = ti.get_previous_ti(state=State.SUCCESS)
        if previous is not None:
            result = previous.xcom_pull()
            if result is not None:
                if result.last_execution_date < start:
                    start = result.last_execution_date


        xcom_q = session.query(XCom).filter(XCom.execution_date >= start).filter(XCom.execution_date < stop)
        xcom_q = xcom_q.filter(XCom.dag_id == self.target_dag_id).filter(XCom.task_id == self.target_task_id)

        xcom_q = xcom_q.order_by(XCom.execution_date.asc(), XCom.timestamp.asc())
        xcom_q = xcom_q.with_entities(XCom.value, XCom.execution_date)

        for xcom in xcom_q.all():
            value = XCom.deserialize_value(xcom)
            if isinstance(value, list):
                if result is None and len(value) > 1:
                    result = RollupState(
                        result=reduce(self.concat, value[1:], value[0]),
                        last_execution_date=xcom.execution_date,
                    )
                else:
                    result = RollupState(
                        result=reduce(self.concat, value, result.result),
                        last_execution_date=xcom.execution_date,
                    )
            else:
                if result is not None:
                    result = RollupState(
                        result=self.concat(value, result.result),
                        last_execution_date=xcom.execution_date,
                    )
                else:
                    result = RollupState(
                        result=value,
                        last_execution_date=xcom.execution_date,
                    )

        return result