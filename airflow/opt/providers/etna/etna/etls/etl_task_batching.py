from datetime import datetime
from typing import Tuple

from airflow.models.taskinstance import Context, TaskInstance
from airflow.utils.timezone import utc

batch_start_context_key = 'prev_data_interval_end_success'
batch_end_context_key = 'data_interval_end'

# A theoretical floor for date ranges.
LOWEST_BOUND = datetime(2010, 1, 1).replace(tzinfo=utc)

def get_batch_range(context: Context) -> Tuple[datetime, datetime]:
    start_date = context[batch_start_context_key] or LOWEST_BOUND
    end_date = context[batch_end_context_key]
    return start_date, end_date