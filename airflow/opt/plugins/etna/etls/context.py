from datetime import datetime, timedelta
from typing import Tuple, Optional

import dateutil.parser
from airflow.models import Variable
from airflow.models.taskinstance import Context, TaskInstance
from airflow.utils.timezone import utc

batch_start_context_key = 'prev_data_interval_end_success'
batch_end_context_key = 'data_interval_end'

# A theoretical floor for date ranges.
LOWEST_BOUND = datetime(2010, 1, 1).replace(tzinfo=utc)

def complete_batch_range(context: Context):
    start, end = get_batch_range(context)
    ti: TaskInstance = context['ti']
    Variable.set(f"{ti.dag_id}-{ti.task_id}-lower-bound", end.isoformat())

def get_batch_range(context: Context) -> Tuple[datetime, datetime]:
    ti: TaskInstance = context['ti']

    # Tasks can be added to a dag after the fact, but need to perform catchup broader than the dag's own lower bound.
    # We keep a running variable of task specific lower bounds in order to manage this.
    task_lower_bound: Optional[str] = Variable.get(f"{ti.dag_id}-{ti.task_id}-lower-bound", None, deserialize_json=False)
    start_date: Optional[datetime] = None
    if task_lower_bound is not None:
        task_lower_bound_dt: datetime = dateutil.parser.isoparse(task_lower_bound)
        start_date = context[batch_start_context_key]
        if start_date is not None and start_date > task_lower_bound_dt:
            start_date = task_lower_bound_dt

    if start_date is None:
        start_date = LOWEST_BOUND


    start_date = start_date.replace(microsecond=0)
    end_date = context[batch_end_context_key]

    # For etls where time boundaries should overlap slightly to ensure precision isn't a factor in missed
    # entries.
    end_date += timedelta(seconds=1)

    end_date = end_date.replace(microsecond=0)

    return start_date, end_date
