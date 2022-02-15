from datetime import datetime
from typing import Tuple

from airflow.models.taskinstance import Context, TaskInstance
from airflow.utils.timezone import utc

batch_start_context_key = 'data_interval_start'
batch_end_context_key = 'data_interval_end'

# A theoretical floor for date ranges.
LOWEST_BOUND = datetime(2010, 1, 1).replace(tzinfo=utc)

def get_batch_range(context: Context) -> Tuple[datetime, datetime]:
    ti: TaskInstance = context['ti']
    if ti.get_previous_dagrun() is None:
        start_date = LOWEST_BOUND
    else:
        start_date = context[batch_start_context_key]
    end_date = context[batch_end_context_key]
    return start_date, end_date