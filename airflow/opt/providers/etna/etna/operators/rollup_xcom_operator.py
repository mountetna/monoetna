from functools import reduce
from typing import Any, Optional, Callable

from airflow.models import BaseOperator, XCom
from airflow.utils.session import provide_session
from etna.etls.etl_task_batching import get_batch_range
from sqlalchemy.orm import Session


class RollupXcomOperator(BaseOperator):
    concat: Callable[[Any, Any], Any]
    target_dag_id: str
    target_task_id: str

    def __init__(self,
                 target_dag_id: str,
                 target_task_id: str,
                 concat: Callable[[Any, Any], Any],
                 task_id: Optional[str] = None, **kwargs
                 ):
        super().__init__(task_id or f"rollup_{target_dag_id}_{target_task_id}", **kwargs)
        self.concat = concat
        self.target_dag_id = target_dag_id
        self.target_task_id = target_task_id

    @provide_session
    def execute(self, context: Any, session: Session = None):
        start, stop = get_batch_range(context)
        xcom_q = session.query(XCom).filter(XCom.execution_date >= start).filter(XCom.execution_date < stop)
        xcom_q = xcom_q.filter(XCom.dag_id == self.target_dag_id).filter(XCom.task_id == self.target_task_id)

        xcom_q = xcom_q.order_by(XCom.execution_date.asc(), XCom.timestamp.asc())
        xcom_q = xcom_q.with_entities(XCom.value)

        result = None

        for xcom in xcom_q.all():
            value = XCom.deserialize_value(xcom)
            if isinstance(value, list):
                if result is None and len(value) > 1:
                    result = reduce(self.concat, value[1:], value[0])
                else:
                    result = reduce(self.concat, value, result)

        return result